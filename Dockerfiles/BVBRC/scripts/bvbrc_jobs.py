import os, subprocess, time, re, json, argparse, logging
from bs4 import BeautifulSoup
from datetime import date


class RunSubprocess:
    """
    Runs a shell command via pythons subprocess module. Stores the completed process instance and the input command.

    Methods:
        get_output(self, length):
            Used to retrieve the output of a command from the completed process instance stdout.

            Args:
                length (optional- int): Length, in bytes, of returned output. Defaults to 1024.

        str_command(self):
            Returns a formatted version of the input command.

            Args: None

        run(self):
            Creates the subprocess and runs the command. This will be called automatically unless the instance was created with auto_run=False
    """

    def __init__(self, command: list[str], auto_run: bool = True, timeout: bool = None):
        self.command = command
        self.timeout = timeout
        self.auto_run = auto_run
        self.process = None

        if self.auto_run:
            self.run()

    # TODO: Add a timeout to this?
    def run(self):
        logging.debug(f"Running command {self.str_command()}")
        self.process = subprocess.run(
            self.command, capture_output=True, text=True, timeout=self.timeout
        )
        logging.debug(f"Command ran with output:\n{self.get_output()}")

    def get_output(self):
        return self.process.stdout

    def str_command(self):
        return " ".join(self.command)

    def __str__(self) -> str:
        # Python magic/dunder method. Called when the object is printed.
        return f"Command: {self.str_command()}\nOutput: {self.get_output()}\nself.timeout: {self.timeout}\nself.auto_run: {self.auto_run}"


def bvbrc_ls_files(dir) -> list[str] | None:
    """
    Runs the p3-ls command on the given BV-BRC workspace path.

    Args:
        dir (string): BV-BRC workspace path

    Returns:
        list[str]: list of files in the directory (if it exists) or None
    """
    ls_command = RunSubprocess(command=["p3-ls", dir, "--one-column"])
    ls_out = ls_command.get_output()
    if "Object not found" in ls_out:
        logging.warning(f"The file or directory {dir} does not exist on BV-BRC.")
        return None

    # splits files into a list()
    # the "--one-column" option uses the newline char to separate the filenames
    # this means that any files with spaces in their names will not break
    files = ls_out.split("\n")
    for file in files:
        file.strip("\r")
    return files


def bvbrc_job_running(job_id: int) -> bool:
    """
    Checks if a BV-BRC job is still running by querying its status.

    Args:
        job_id (int): BV-BRC job ID.

    Raises:
        TypeError: If job_id is None or not an integer.
        ValueError: If the job ID is invalid or not found.

    Returns:
        bool: True if the job is still running, otherwise False.
    """
    if job_id is None:
        # Log a warning if no job ID was provided, indicating a possible submission issue.
        logging.warning("The job id is None. Check that job was properly submitted")
        raise TypeError # Raise an error as None is not a valid job ID

    # Run the BV-BRC job status command to check the current state of the job
    job_status = RunSubprocess(command=["p3-job-status", job_id]).get_output()
    print(job_status) # Print the job status for debugging purposes

    if "job not found" in job_status:
        # Log an issue if the job ID is invalid or was not recorded properly.
        logging.debug(
            f"Invalid Job Id of {job_id}: Check the get_job_id() func for valid output."
        )
        raise ValueError # Raise an error for an invalid job ID
    
    # Determine if the job is still running by checking for "in-progress" status.
    # return "in-progress" in job_status
    return any(status in job_status for status in {"in-progress", "queued", "pending"})


def get_job_id(command_output: str) -> int:
    """
    Checks the output of a Bv-BRC "submit job" command to see if there was a
    job id returned.

    Finds id by checking for string "id " in return.

    Args:
        command_output (str): String output of a submit job command (e.g. p3-submit-genome-assembly)

    Returns:
        int: BV-BRC job id
    """
    match = re.search(r"id (\d+)", command_output)
    if match:
        job_id = match.group(1)
        return job_id
    return None


def genome_assembly_job(
    read1: str | os.PathLike, read2: str | os.PathLike, dry_run: bool = False
) -> None:
    """
    Must be logged into the BV-BRC CLI.
    Runs a genome assembly job using BV-BRC.

    This function:
    1. Creates the required directory structure and uploads input files.
    2. Submits a genome assembly job (optionally as a dry run).
    3. Processes and extracts results from the job output.

    Args:
        read1 (str | os.PathLike): Local path to read1.
        read2 (str | os.PathLike): Local path to read2.
        dry_run (bool, optional): If True, runs the job without actual execution. Defaults to False.
    """

    def __generate_directory_structure() -> tuple[str, str]:
        """
        Creates the BV-BRC workspace directory structure, uploads input files, 
        and initializes an output directory.

        Returns:
            tuple[str, str]: Paths to the uploaded read1 and read2 files in the workspace.
        """

        directory_commands = [
            ["p3-mkdir", RAW_READS_DIR],
            ["p3-mkdir", ASSEMBLY_DIR],
            ["p3-cp", read1, f"ws:{RAW_READS_DIR}"],
            ["p3-cp", read2, f"ws:{RAW_READS_DIR}"],
            ["mkdir", "bvbrc_asm_output"],
        ]

        for command in directory_commands:
            RunSubprocess(command) # Execute each command

        print(RunSubprocess(["ls", "bvbrc_asm_output"]).get_output()) # List output directory contents
        
        raw_reads_dir_ls = bvbrc_ls_files(RAW_READS_DIR) # Retrieve uploaded filenames
        read1_ws_name, read2_ws_name = raw_reads_dir_ls[0], raw_reads_dir_ls[1] # Get workspace filenames
        
        logging.debug("Created directory structure")

        return (
            f"ws:{RAW_READS_DIR}/{read1_ws_name}",
            f"ws:{RAW_READS_DIR}/{read2_ws_name}",
        )

    def bvbrc_genome_assembly(ws_read1, ws_read2, dry_run: bool) -> None:
        """
        Submits a genome assembly job to BV-BRC.

        Args:
            ws_read1 (str): Path to read1 in the workspace.
            ws_read2 (str): Path to read2 in the workspace.
            dry_run (bool): If True, performs a dry run without execution.
        """
        if dry_run:
            logging.warning("Running dry run job. No data will be submitted.")
            assembly_job = RunSubprocess(
                command=[
                    "p3-submit-genome-assembly",
                    "--trim-reads",
                    "--workspace-upload-path",
                    RAW_READS_DIR,
                    "--recipe unicycler",
                    "--paired-end-lib",
                    "--dry-run",
                    ws_read1,
                    ws_read2,
                    "--min-contig-len",
                    "300",
                    "--min-contig-cov",
                    "30",
                    f"ws:{ASSEMBLY_DIR}",
                    args.sample_name,
                ]
            )
            return
        # Submit real genome assembly job
        assembly_job = RunSubprocess(
            command=[
                "p3-submit-genome-assembly",
                "--trim-reads",
                "--workspace-upload-path",
                RAW_READS_DIR,
                "--recipe unicycler",
                "--paired-end-lib",
                ws_read1,
                ws_read2,
                "--min-contig-len",
                "300",
                "--min-contig-cov",
                "30",
                f"ws:{ASSEMBLY_DIR}",
                args.sample_name,
            ]
        )

        job_id = get_job_id(assembly_job.get_output()) # Extract job ID from output
        if job_id is not None:
            while bvbrc_job_running(job_id): # Monitor job status
                time.sleep(60)  # Wait before checking again

    def __handle_asm_output() -> None:
        """
        Processes the genome assembly output and extracts relevant metrics.
        """
        output_commands = [
            ["touch", "bvbrc_asm_output/output_path.txt"], # Create an output file
        ]

        for command in output_commands:
            RunSubprocess(command) # Execute each command

        for file in ASM_OUTPUTS.values():
            command = RunSubprocess(["p3-cp", file, "bvbrc_asm_output"]) # Copy output files

        # Parse the assembly report
        with open(
            f"bvbrc_asm_output/{args.sample_name}_AssemblyReport.html", "r"
        ) as file:
            soup = BeautifulSoup(file, "html.parser")

        # Extract preprocessing metrics - Find the table with the header "Preprocessing of Reads"
        table = soup.find("h2", string="Preprocessing of Reads").find_next("table")

        # Extract the number of reads from the table
        table_data = []
        num_reads = None
        rows = table.find_all("tr")

        for i, row in enumerate(rows):
            if i == 0:
                cells = row.find_all("th")
            else:
                cells = row.find_all("td")

            cells_data = list(map(lambda x: x.get_text(), cells))
            table_data.append(cells_data)

        num_reads = table_data[1][table_data[0].index("Num Reads")]
        average_size = table_data[1][table_data[0].index("Avg Length")]

        # Extract assembly metrics - Find the table with the header "Assembly"
        table2 = soup.find("h2", string="Assembly").find_next("table")
        rows2 = table2.find("tbody").find_all("tr")
        fasta_file_size = None
        for row in rows2:
            cell = row.find_all("td")
            if cell[0].get_text() == "contigs.fasta file size:":
                fasta_file_size = cell[1].get_text()

        # Extract additional details from JSON output
        with open(f"bvbrc_asm_output/{args.sample_name}_run_details.json", "r") as f:
            data = json.load(f)
            average_depth = data["contig_filtering"]["average depth (short reads)"]
            contigs_above_threshold = data["contig_filtering"][
                "num contigs above thresholds"
            ]
            contigs_below_threshold = data["contig_filtering"][
                "num contigs below thresholds"
            ]

        # Write extracted details to output file
        with open("bvbrc_asm_output/output_path.txt", "w") as f:
            f.write(f"Contigs Workspace Path: {ASM_OUTPUTS['CONTIG_OUTPUT']}\n")
            f.write(f"Read Timestamp: {FROM_EPOCH}\n")
            f.write(f"Contig.fasta File Size: {fasta_file_size}\n")
            f.write(f"Number of Reads: {num_reads}\n")
            f.write(f"Average Read Length: {average_size}\n")
            f.write(f"Average Read Depth: {average_depth}\n")
            f.write(f"Number of Contigs Above Threshold: {contigs_above_threshold}\n")
            f.write(f"Number of Contigs Below Threshold: {contigs_below_threshold}\n")

    # Execute the workflow steps
    ws_read1, ws_read2 = __generate_directory_structure() # Step 1: Prepare workspace
    bvbrc_genome_assembly(ws_read1, ws_read2, dry_run) # Step 2: Run genome assembly
    __handle_asm_output() # Step 3: Process output


def bvbrc_annotation_analysis(
    contigs: str | os.PathLike,
    output_path: str | os.PathLike,
    scientific_name: str,
    sample_name: str,
    taxonomy_id: str = "2",
) -> None:
    """
    Submits assembly contig to BVBRC for annotation and assembly.

    Args:
        sample_name (str): Name of sample for analysis
        contigs (str | os.PathLike): assembly file
        output_path (str | os.PathLike): BVBRC path to upload results to
        output_name (str): Name of output file
        scientific_name (str): Scientific name of sample
        taxonomy_id (int, optional): Defaults to 2.
    """

    # Define the output name using the sample name
    output_name = f"{sample_name}_output"

    # Create the output directory in BV-BRC workspace if it doesn't already exist
    # NOTE: This is a workspace path, but should not be prefixed with "ws:"
    RunSubprocess(["p3-mkdir", output_path])

    # Construct the command to submit the contigs for annotation via BV-BRC
    cga_command = [
        "p3-submit-CGA",
        "-contigs",
        contigs,
        "-scientific-name",
        scientific_name,
        "-taxonomy-id",
        str(taxonomy_id),
        "-code",
        "11",
        "-domain",
        "Bacteria",
        output_path,
        output_name,
    ]

    # Run the annotation submission command
    cga_job = RunSubprocess(command=cga_command)

    # Extract the job ID from the command output
    job_id = get_job_id(cga_job.get_output())

    # Wait for the annotation job to complete, checking its status periodically
    count = 0 
    while bvbrc_job_running(job_id):
        print(count)
        count += 1
        logging.debug(f"Waiting for BV-BRC annotation job to finish. Attempt {count}")
        time.sleep(60)  # This timing may need to be adjusted.

    # Issue: If the expected output files do not exist or change location,
    #        the command may hang indefinitely.
    # Solution: Consider adding a timeout or error handling for missing files.

    # Commands to retrieve annotation output files from BV-BRC workspace
    output_commands = [
        ["mkdir", "bvbrc_cga_output"],
        ["p3-cp", CGA_OUTPUTS["FULL_GENOME_REPORT"], "bvbrc_cga_output"],
        ["p3-cp", CGA_OUTPUTS["ANNOTATION_GENOME"], "bvbrc_cga_output"],
        ["p3-cp", CGA_OUTPUTS["ANNOTATION_QUALITY"], "bvbrc_cga_output"],
        ["p3-cp", CGA_OUTPUTS["ANNOTATION_RESISTANCE"], "bvbrc_cga_output"],
        ["p3-cp", CGA_OUTPUTS["ANNOTATION_PROTEIN"], "bvbrc_cga_output"],
    ]

    # Execute each command to retrieve the annotation results
    for command in output_commands:
        logging.debug(f"Running Command: {command}")
        RunSubprocess(command)


def bvbrc_annotation_analysis_local(
    contigs_local_path, output_path, scientific_name, sample_name, taxonomy_id=2
):
    """
    Runs a BV-BRC annotation analysis locally by first uploading contigs to the workspace,
    and then calling the annotation function.

    Args:
        contigs_local_path (str): Local path to the contigs file.
        output_path (str): Directory where the annotation results should be saved.
        scientific_name (str): Scientific name of the sample.
        sample_name (str): Name of the sample.
        taxonomy_id (int, optional): Taxonomy ID of the sample. Defaults to 2.

    Raises:
        ValueError: If no assembly files are found in the workspace after upload.
    """
    # Commands to create a workspace directory and upload contigs to BV-BRC workspace
    workspace_directory_commands = [
        ["p3-mkdir", ASSEMBLY_DIR],
        ["p3-cp", contigs_local_path, f"ws:{ASSEMBLY_DIR}"],
    ]
    # Execute each command to set up the workspace
    for command in workspace_directory_commands:
        RunSubprocess(command)

    # List files in the assembly directory to verify upload
    assembly_files = bvbrc_ls_files(ASSEMBLY_DIR)

    if assembly_files is None:
        # Raise an error if no files are found after upload
        raise ValueError

    # Construct the workspace path to the uploaded contigs file
    contigs_ws_path = f"ws:{ASSEMBLY_DIR}/{assembly_files[0]}"

    # Run the BV-BRC annotation analysis using the uploaded contigs
    bvbrc_annotation_analysis(
        contigs=contigs_ws_path,
        output_path=output_path,
        scientific_name=scientific_name,
        sample_name=sample_name,
        taxonomy_id=taxonomy_id,
    )


parser = argparse.ArgumentParser(
    prog="AutoBVBRC",
    description="What the program does",
    epilog="Text at the bottom of help",
)

# JOB SELECTION
selected_job = parser.add_mutually_exclusive_group()
selected_job.add_argument("-asm", "--genome-assembly", action="store_true")
selected_job.add_argument("-cga", "--complete-genome-analysis", action="store_true")
selected_job.add_argument(
    "-cgal", "--complete-genome-analysis-local", action="store_true"
)
selected_job.add_argument("-dev", "--development", action="store_true")

# Logging
command_output = parser.add_mutually_exclusive_group()
command_output.add_argument(
    "-q",
    "--quiet",
    action="store_true",
    help="Skips all logging to stdout except CRITICAL errors",
)
command_output.add_argument(
    "-v", "--verbose", action="store_true", help="More comprehensive logging to stdout"
)
command_output.add_argument(
    "-d",
    "--debug",
    action="store_true",
    help="Most comprehensive logging to stdout for DEBUGGING",
)

# COMMON
parser.add_argument(
    "-u",
    "--username",
    type=str,
    metavar="",
    required=False,
    help="BV-BRC username for account login",
)
parser.add_argument(
    "-p",
    "--password",
    type=str,
    metavar="",
    required=False,
    help="BV-BRC password for account login",
)
parser.add_argument(
    "-n",
    "--sample-name",
    type=str,
    metavar="",
    required=False,
    help="Name of bacterial sample",
)
parser.add_argument(
    "--dry-run",
    action="store_true",
    help="A dry run will validate input, but not submit a BVBRC job",
)

# CGA
parser.add_argument(
    "-a",
    "--assembly-file",
    type=str,
    metavar="",
    required=False,
    help="BV-BRC workspace path of assembly",
)
parser.add_argument(
    "-t",
    "--timestamp",
    type=str,
    metavar="",
    required=False,
    help="Timestamp of the reference assembly job",
)
parser.add_argument(
    "-sci",
    "--scientific-name",
    type=str,
    metavar="",
    required=False,
    help="Scientific name of bacterial sample",
)
parser.add_argument(
    "-tax",
    "--taxonomy-id",
    type=int,
    metavar="",
    required=False,
    default=2,
    help="Taxonomy ID of sample, defaults to 2",
)

# ASM
parser.add_argument(
    "-r1",
    "--read1",
    type=str,
    metavar="",
    required=False,
    default=2,
    help="Read 1 file for assembly",
)
parser.add_argument(
    "-r2",
    "--read2",
    type=str,
    metavar="",
    required=False,
    default=2,
    help="Read 2 file for assembly",
)

args = parser.parse_args()


class Log:
    HELP_TIP = f"Run Command with -h or --help to see a complete list of options."


if args.quiet:
    log_level = logging.CRITICAL
elif args.verbose:
    log_level = logging.INFO
elif args.debug:
    log_level = logging.DEBUG
else:
    log_level = logging.WARNING


logging.basicConfig(
    level=log_level,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

########## CONSTANTS ##########
FROM_EPOCH = int(time.time())

PROJECT_ROOT = f"/{args.username}@bvbrc/home/CAMRA_WDL"

DATE_DIR = f"{PROJECT_ROOT}/{date.today()}"
SAMPLE_DIR = f"{DATE_DIR}/{FROM_EPOCH}_{args.sample_name}"
RAW_READS_DIR = f"{SAMPLE_DIR}/raw_reads"
ASSEMBLY_DIR = f"{SAMPLE_DIR}/assembly"

ASM_OUTPUT_DIR = f"ws:{ASSEMBLY_DIR}/.{args.sample_name}"
ASM_DETAILS_DIR = f"ws:{ASSEMBLY_DIR}/.{args.sample_name}/details"

# ASM OUTPUTS
ASM_OUTPUTS = {
    "CONTIG_OUTPUT": f"{ASM_OUTPUT_DIR}/{args.sample_name}_contigs.fasta",
    "ASSEMBLY_REPORT": f"{ASM_OUTPUT_DIR}/{args.sample_name}_AssemblyReport.html",
    "BANDAGE_PLOT": f"{ASM_DETAILS_DIR}/{args.sample_name}_assembly_graph.plot.svg",
    "DETAILS_JSON": f"{ASM_DETAILS_DIR}/{args.sample_name}_run_details.json",
}

cga_output_path = f"{DATE_DIR}/{args.timestamp}_{args.sample_name}/CGA"

CGA_OUTPUTS = {
    "FULL_GENOME_REPORT": f"ws:{cga_output_path}/.{args.sample_name}_output/FullGenomeReport.html",
    "ANNOTATION_GENOME_REPORT": f"ws:{cga_output_path}/.{args.sample_name}_output/.annotation/GenomeReport.html",
    "ANNOTATION_XLS": f"ws:{cga_output_path}/.{args.sample_name}_output/.annotation/annotation.xls",
    "ANNOTATION_GENOME": f"ws:{cga_output_path}/.{args.sample_name}_output/.annotation/quality.json",
    "ANNOTATION_QUALITY": f"ws:{cga_output_path}/.{args.sample_name}_output/.annotation/annotation.genome",
    "ANNOTATION_RESISTANCE": f"ws:{cga_output_path}/.{args.sample_name}_output/.annotation/load_files/genome_amr.json",
    "ANNOTATION_PROTEIN": f"ws:{cga_output_path}/.{args.sample_name}_output/.annotation/annotation.feature_protein.fasta",
}

########## END CONSTANTS ##########


def main(args):
    """
    Main function to execute different genomic analysis workflows based on user-provided arguments.

    Args:
        args: Parsed command-line arguments containing workflow options and input parameters.
    """

    # Warn the user if they are running the script in dry-run mode
    if args.dry_run:
        logging.warning(
            "You submitted a job with the '--dry-run' option. This may break the output of the workflow if the output has not already been generated."
        )

    # Run genome assembly if the corresponding argument is provided
    if args.genome_assembly:
        genome_assembly_job(read1=args.read1, read2=args.read2, dry_run=args.dry_run)

    # Run complete genome analysis using the BV-BRC annotation pipeline
    elif args.complete_genome_analysis:
        bvbrc_annotation_analysis(
            contigs=args.assembly_file,
            output_path=cga_output_path,
            scientific_name=args.scientific_name,
            sample_name=args.sample_name,
            taxonomy_id=args.taxonomy_id,
        )

    # Run local complete genome analysis using a local path for input contigs
    elif args.complete_genome_analysis_local:
        bvbrc_annotation_analysis_local(
            contigs_local_path=args.assembly_file,
            output_path=cga_output_path,
            scientific_name=args.scientific_name,
            sample_name=args.sample_name,
            taxonomy_id=args.taxonomy_id,
        )

    # Development mode: Logs basic messages for testing purposes
    elif args.development:
        logging.warning("Dev Run Started")
        logging.warning("Dev Run Finished")

    # If no valid job type is provided, log an error message
    else:
        logging.error("Error in the job selection. Function - main()")

# Ensure the script runs only when executed directly (not when imported)
if __name__ == "__main__":
    main(args)
