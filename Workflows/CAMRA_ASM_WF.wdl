version 1.0

import "../Tasks/BV-BRC_tasks.wdl" as bvbrc

workflow assembly_workflow {
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
    call bvbrc.run_genome_assembly {
        input:
            read1 = read1,
            read2 = read2,
            sample_name = sample_name,
            username = BVBRC_user,
            password = BVBRC_password,

    }

    output {
        String testOut = run_genome_assembly.testOut
    }

}