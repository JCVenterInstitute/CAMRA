version 1.0

task run_fastQC {
    meta {
        # TODO
    }

    input {
        File read1
        File read2
    }

    runtime {
        docker: 'staphb/fastqc:0.12.1'
    }

    command <<<
        echo "fastQC started"
        mkdir fastQC_output
        fastqc -o fastQC_output -j /usr/bin/java -f fastq ~{read1} ~{read2}
        echo "fastQC ended"
    >>>

    output {
        # TODO
        # date
        # version
        # File fastQC_report = "fastQC_output/report.txt"
    }
}