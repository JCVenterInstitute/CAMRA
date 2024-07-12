import os, io, sys, pty, subprocess, time, re

class PseudoTerminalPair():
    """
    Creates pseudo terminal pair with a command to run.

    Methods:
        run_process(self, command):
            Creates a subprocess in the worker terminal. Runs given command and waits for completion.

            Args:
                command (list[str]): Bash command broken into args as a list of strings

                
        get_output(self, length):
            Used to retrieve the output of a command. Decoded to a string.

            Args:
                length (optional- int): Length, in bytes, of returned output. Defaults to 1024.

    """
    
    def __init__(self, command: list[str]):
        self.controller, self.worker = pty.openpty()
        self.run_process(command)
    
    def run_process(self, command):
        try:
            process = subprocess.Popen(command, stdin=self.worker, 
                            stdout=self.worker, stderr=self.worker, 
                            shell=False)
            process.wait()
        except Exception as e:
            print(f"There was an error running the command '{command}':\n{e}")
    
    def get_controller(self):
        return self.controller
    
    def get_worker(self):
        return self.worker
    
    def get_output(self, length = 1024) -> str:
        output = os.read(self.controller, length).decode()
        return output
    
    def close_terminal(self):
        os.close(self.controller)
        os.close(self.worker)

def run_command(command:list[str], print_out:bool = False) -> None:
    """
    Generates a pseudo terminal pair and executes bash command.
    Automatically kills terminal upon completion.

    Args:
        command (list[str]): Bash command broken into args as a list of strings
        print_out (optional- bool): If True, prints command output to stdout.
    """
    run_command = PseudoTerminalPair(command=command)
    if print_out:
        print(run_command.get_output())
    run_command.close_terminal()


def bvbrc_job_running(job_id) -> bool:
    if job_id is None:
        raise TypeError
    
    job_status_process = PseudoTerminalPair(command=['p3-job-status', job_id])
    status = 'in-progress' in job_status_process.get_output()
    job_status_process.close_terminal()
    return status

def get_job_id(command_output:str) -> int:
    # NOTE:
    # This is a slightly brittle method. It relies on the output of the command including
    # the job id prefixed by "id". Would break if bvbrc changed the output of their 
    # assembly command.
    match = re.search(r'id (\d+)', command_output)
    if match:
        job_id = match.group(1)
        return job_id
    return None



def bvbrcGenomeAssembly(user: str, sample_name: str, read1: str | os.PathLike, read2: str | os.PathLike ) -> None:
    """
    Submits two read files to BV-BRC for assembly. 
    Must be logged into the BV-BRC cli.

    Args:
        read1 (str | os.PathLike): Path to read1
        read2 (str | os.PathLike): Path to read2
    """
    root_dir = f"/{user}@bvbrc/home"
    sample_dir = f"{root_dir}/assemblies/{sample_name}"
    raw_reads_dir = f"{sample_dir}/raw_reads"
    assembly_dir = f"{sample_dir}/assembly"
    output_path = f"ws:{assembly_dir}/.{sample_name}/{sample_name}_contigs.fasta"
    
    ls_command = PseudoTerminalPair(command=['p3-ls', sample_dir])
    if 'Object not found' in ls_command.get_output():
        directory_commands = [['p3-mkdir', raw_reads_dir],
                          ['p3-mkdir', assembly_dir],
                          ['p3-cp', read1, f"ws:{raw_reads_dir}"],
                          ['p3-cp', read2, f"ws:{raw_reads_dir}"],
                          ['mkdir', 'bvbrc_asm_output'],]
        for command in directory_commands:
            run_command(command)
        
        # NOTE:
        # Reading from the ls command is generally considered to be bad practice.
        # This could cause errors if there are unusual file names (filepaths with newline characters, etc.).
        # BVBRC doesn't seem to have a glob method or some better way of doing this.
        # However, it should work fine in the context of this project as we define filepaths internally.
        raw_ls = PseudoTerminalPair(command=['p3-ls', f"{raw_reads_dir}/", '--one-column'])
        raw_ls_output = raw_ls.get_output()
        raw_files = raw_ls_output.split('\n')
        read1_name, read2_name = raw_files[0].strip('\r'), raw_files[1].strip('\r')
            
        assembly_job = PseudoTerminalPair(command=['p3-submit-genome-assembly', '--trim-reads', 
                            '--workspace-upload-path', raw_reads_dir, 
                            '--recipe unicycler','--paired-end-lib', 
                            f"ws:{raw_reads_dir}/{read1_name}", f"ws:{raw_reads_dir}/{read2_name}",  '--min-contig-len', 
                            '300', '--min-contig-cov', '30', f"ws:{assembly_dir}", sample_name])
        
        # Get the job id
        job_id = get_job_id(assembly_job.get_output())
        
        while bvbrc_job_running(job_id):
            time.sleep(10)  # This timing may need to be adjusted.
        
        print(f"Job {job_id} Finished")
        
        
        # Close all of the pseudo terminals
        raw_ls.close_terminal()
        assembly_job.close_terminal()
       
    else:
        print("Directory already exists:")

    output_commands = [['mkdir', 'bvbrc_asm_output'],
                        ['touch', 'bvbrc_asm_output/output_path.txt'],
                        ['p3-cp', output_path, 'bvbrc_asm_output']]
    
    for command in output_commands:
        run_command(command)

    with open('bvbrc_asm_output/output_path.txt', 'w') as f:
        print(output_path, file=f)


def bvbrcGenomeAnnotation(user: str, sample_name: str, assembly_filepath: str | os.PathLike, scientific_name: str=None) -> None:#

    # TODO change asm dir to user input
    root_dir = f"/{user}@bvbrc/home"
    sample_dir = f"{root_dir}/assemblies/{sample_name}"
    assembly_dir = f"{sample_dir}/assembly"
    output_dir = f"{sample_dir}/annotation"

    run_command(['p3-mkdir', output_dir])
    
    # TODO allow for user to provide file for upload or bvbrc location of file
    annotation_command = ['p3-submit-genome-annotation', '--workspace-path-prefix', assembly_dir, '-f',
                        '--contigs-file', assembly_filepath, '-d', 'Bacteria']
    if scientific_name:
        annotation_command.append('-n')
        annotation_command.append(scientific_name)
    
    annotation_command.append(output_dir)
    annotation_command.append(f"{sample_name}_annotation")

    annotation_job = PseudoTerminalPair(command=annotation_command)

    job_id = get_job_id(annotation_job.get_output())

    while bvbrc_job_running(job_id):
        time.sleep(10)  # This timing may need to be adjusted.
    
    print(f"Job {job_id} Finished")

    annotation_job.close_terminal()

def completeGenomeAnalysis(user: str, sample_name: str, assembly_filepath: str | os.PathLike, scientific_name: str=None):
    pass
    
def main():
    match sys.argv[1].lower():

        case 'assembly' | 'asm':
            bvbrcGenomeAssembly(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

        case 'annotation' | 'annot':
            match len(sys.argv):
                case 5:
                    bvbrcGenomeAnnotation(sys.argv[2], sys.argv[3], sys.argv[4])
                case 6:
                    bvbrcGenomeAnnotation(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
                case _:
                    raise IndexError
                
        case 'analysis' | 'cga':
            pass

        case _:
            raise ValueError




if __name__ == "__main__":
    main()