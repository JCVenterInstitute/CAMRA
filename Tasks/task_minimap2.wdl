version 1.0

task run_minimap2 {
    input{
    File fastq_file
    File reference_fasta
    }

    command <<<

        minimap2 --version | tee VERSION

        minimap2 \
        -ax map-ont \
        ~{reference_fasta} \
        ~{fastq_file} \
        -o minimap2_output.sam

        samtools sort \
        minimap2_output.sam \
        -o minimap2_sorted.sam 
    >>>

    output {
        File minimap2_sorted_output = "minimap2_sorted.sam"
        String minimap2_version = read_string("VERSION")
        
    }

    runtime {
        docker: "nanozoo/minimap2:2.28--9e3bd01"
    }
}