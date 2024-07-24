version 1.0
# All taks relating to BV-BRC
task run_genome_assembly {
    meta {
    author: "Andrew LaPointe"
    email: "andrewrlapointe@gmail.com"
    description: "Task for running the BV-BRC genome assembly tool."
    version: "1.0"
    dockerhub: "https://hub.docker.com/repository/docker/andrewrlapointe/bvbrc/general"
    }

    input {
        File read1
        File read2
        String sample_name
        String username
        String password
    }

    runtime {
        docker: 'andrewrlapointe/bvbrc:3.2'

    }

    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py assembly ~{username} ~{sample_name} ~{read1} ~{read2}
                
        # Extract values
        contigs_workspace_path=$(grep -oP '(?<=Contigs Workspace Path: ).*' bvbrc_asm_output/output_path.txt)
        num_reads=$(grep -oP '(?<=Number of Reads: ).*' bvbrc_asm_output/output_path.txt)
        average_read_depth=$(grep -oP '(?<=Average Read Depth: ).*' bvbrc_asm_output/output_path.txt)
        contigs_above_threshold=$(grep -oP '(?<=Number of Contigs Above Threshold: ).*' bvbrc_asm_output/output_path.txt)
        contigs_below_threshold=$(grep -oP '(?<=Number of Contigs Below Threshold: ).*' bvbrc_asm_output/output_path.txt)
        average_read_length=$(grep -oP '(?<=Average Read Length: ).*' bvbrc_asm_output/output_path.txt)
        contig_fasta_file_size=$(grep -oP '(?<=Contig.fasta File Size: ).*' bvbrc_asm_output/output_path.txt)

        # Save the variables to separate output files
        echo "$contigs_workspace_path" > contigs_workspace_path.txt
        echo "$num_reads" > num_reads.txt
        echo "$average_read_depth" > average_read_depth.txt
        echo "$contigs_above_threshold" > contigs_above_threshold.txt
        echo "$contigs_below_threshold" > contigs_below_threshold.txt
        echo "$average_read_length" > average_read_length.txt
        echo "$contig_fasta_file_size" > contig_fasta_file_size.txt

        rm bvbrc_asm_output/output_path.txt
    >>>

    output {
        File asm_bandage_plot = "bvbrc_asm_output/~{sample_name}_assembly_graph.plot.svg"
        File assembly_file = "bvbrc_asm_output/~{sample_name}_AssemblyReport.html"
        String contigs_workspace_path = read_string("contigs_workspace_path.txt")
        Int contig_fasta_file_size = read_int("contig_fasta_file_size.txt")
        Float average_read_depth = read_float("average_read_depth.txt")
        Int number_reads = read_int("num_reads.txt")
        Float average_read_length = read_float("average_read_length.txt")
        Int contigs_above_threshold = read_int("contigs_above_threshold.txt")
        Int contigs_below_threshold = read_int("contigs_below_threshold.txt")
    }
}


task run_annotation_analysis {
    meta {
        author: "Andrew LaPointe"
        email: "andrewrlapointe@gmail.com"
        description: "Task for running the BV-BRC complete genome analysis tool."
        version: "1.2"
    }

    input {
        String contigs_file
        String username
        String password
        String sample_name
        String scientific_name
        String output_path
        String taxonomy_id
    }

    runtime {
        docker: 'andrewrlapointe/bvbrc:3.3'
    }

    # output path could be changed to be relative to the contigs file location to reduce the number of inputs
    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py cga ~{contigs_file} ~{output_path} ~{sample_name}_output ~{scientific_name} ~{sample_name} ~{taxonomy_id}
    >>>

    output {
        File full_genome_report = "bvbrc_cga_output/FullGenomeReport.html"
        File annotation_genome_report = "bvbrc_cga_output/GenomeReport.html"
        File annotation_xls = "bvbrc_cga_output/annotation.xls"

    }
}