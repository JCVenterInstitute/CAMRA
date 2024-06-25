version 1.0

import "../Tasks/BV-BRC_tasks.wdl" as bvbrc

workflow assembly_workflow {
    meta {
    author: ""
    email: ""
    description: "Create genome assemblyu"
    version: "1.0"
    }

    parameter_meta {
        # TODO
    }

    input {
        File read1
        File read2
        String BVBRC_user
        String BVBRC_password
    }

    # TASKS
    call bvbrc.run_genome_assembly {
        input:
            read1 = read1,
            read2 = read2,
            user_name = BVBRC_user,
            password = BVBRC_password,

    }

    output {
        String testOut = run_genome_assembly.testOut
    }

}