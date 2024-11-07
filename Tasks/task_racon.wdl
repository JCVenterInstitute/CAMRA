version 1.0

task run_racon {
    input{
    File minimap2_sorted_sam
    File nanofilt_fastq
    File reference_fasta
    }

    command <<<
        racon --version | tee VERSION

        racon \
        ~{nanofilt_fastq} \
        ~{minimap2_sorted_sam} \
        ~{reference_fasta} \
        > consensus.fasta 


    >>>

    output {
        File racon_consensus_output = "consensus.fasta"
        String racon_version = read_string("VERSION")
        
    }

    runtime {
        docker: "nanozoo/racon:1.5.0--c32e531"
    }
}