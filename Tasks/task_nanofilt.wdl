version 1.0

task run_nanofilt {
    input{
    File porechop_reads
    }

    command <<<
        NanoFilt --version | tee VERSION
        gunzip -c ~{porechop_reads} | NanoFilt --length 400 --maxlength 8600 -q 15 | gzip > nanofilt.fastq.gz

    >>>

    output {
        File nanofilt_fastq_output = "nanofilt.fastq.gz"
        String nanofilt_version = read_string("VERSION")
    }

    runtime {
        docker: "quay.io/biocontainers/nanofilt:2.7.1--py_0"
    }
}