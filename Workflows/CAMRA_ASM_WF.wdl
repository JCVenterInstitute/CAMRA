version 1.0

import "../Tasks/BV-BRC_tasks.wdl" as bvbrc

workflow genome_assembly {
    meta {
    author: "Andrew LaPointe"
    email: "andrewrlapointe@gmail.com"
    description: "Create genome assembly"
    version: "1.0"
    }

    parameter_meta {
        sample_name     :   "Name of sample isolate"
        read1           :   "raw read 1 fastq.gz or fastq file"
        read2           :   "raw read 2 fastq.gz or fastq file"
    }

    input {
        File read1
        File read2
        String sample_name
        String BVBRC_user
        String BVBRC_password
    }

    # TASKS
    call bvbrc.run_BVBRC_genome_assembly {
        input:
            read1 = read1,
            read2 = read2,
            sample_name = sample_name,
            username = BVBRC_user,
            password = BVBRC_password,
    }

    output {
        File assembly_file = run_BVBRC_genome_assembly.assembly_file
        File asm_bandage_plot = run_BVBRC_genome_assembly.asm_bandage_plot
        String contigs_workspace_path = run_BVBRC_genome_assembly.contigs_workspace_path
        Int contig_fasta_file_size = run_BVBRC_genome_assembly.contig_fasta_file_size
        Int number_reads = run_BVBRC_genome_assembly.number_reads
        Int timestamp = run_BVBRC_genome_assembly.timestamp
        Float average_read_length = run_BVBRC_genome_assembly.average_read_length 
        Float average_read_depth = run_BVBRC_genome_assembly.average_read_depth
        Int contigs_above_threshold = run_BVBRC_genome_assembly.contigs_above_threshold
        Int contigs_below_threshold = run_BVBRC_genome_assembly.contigs_below_threshold
    }

}