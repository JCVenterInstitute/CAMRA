import os, io, sys, pty, subprocess, time, re

# TODO - Make pseudo terminals self closing unless specified.

class PseudoTerminalPair():
    
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
    
    def get_output(self):
        output = os.read(self.controller, 1024).decode()
        return output
    
    def close_terminal(self):
        os.close(self.controller)
        os.close(self.worker)


def bvbrcGenomeAssembly(user: str, sample_name: str, read1: str | os.PathLike, read2: str | os.PathLike ) -> None:
    """_summary_

    Args:
        read1 (str | os.PathLike): Path to read1
        read2 (str | os.PathLike): Path to read2
    """
    root_dir = f"/{user}@bvbrc/home"
    sample_dir = f"{root_dir}/assemblies/{sample_name}"
    raw_reads_dir = f"{sample_dir}/raw_reads"
    assembly_dir = f"{sample_dir}/assembly"    
    
    ls_command = PseudoTerminalPair(command=['p3-ls', sample_dir])
    if 'Object not found' in ls_command.get_output():
        directory_commands = [['p3-mkdir', raw_reads_dir],
                          ['p3-mkdir', assembly_dir],
                          ['p3-cp', read1, f"ws:{raw_reads_dir}"],
                          ['p3-cp', read2, f"ws:{raw_reads_dir}"],
                          ['mkdir', 'bvbrc_asm_output'],]
        for command in directory_commands:
            terminal_obj_list = [PseudoTerminalPair(command) for command in directory_commands]
        
        # NOTE:
        # Reading from the ls command us considered to be bad practice.
        # This could cause errors if there are unusual file names.
        # BVBRC doesn't seem to have a glob method or some better way of doing this.
        # However, it should work fine in the context of this project.
        # Here's to hoping. Cheers
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
        match = re.search(r'id (\d+)', assembly_job.get_output())
        if match:
            job_id = match.group(1)
            
        
        while bvbrc_job_running(job_id):
            time.sleep(10)  # This timing may need to be adjusted.
        
        print(f"Job {job_id} Finished")
        
        
        # Close all of the pseudo terminals
        for i in terminal_obj_list:
            i.close_terminal()
        raw_ls.close_terminal()
        assembly_job.close_terminal()
        
        PseudoTerminalPair(command=['p3-cp', f"{assembly_dir}/{sample_name}", 'bvbrc_asm_output'])
        return
       
    else:
        print("Directory already exists:")
        return


def bvbrc_job_running(job_id):
    
    job_status_process = PseudoTerminalPair(command=['p3-job-status', job_id])
    status = 'in-progress' in job_status_process.get_output()
    job_status_process.close_terminal()
    return status


bvbrcGenomeAssembly(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
