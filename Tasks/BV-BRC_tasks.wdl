version 1.0
# All taks relating to BV-BRC
task run_genome_assembly {
    meta {
    author: ""
    email: ""
    description: "Task for running the BV-BRC genome assembly tool."
    version: "1.0"
    dockerhub: "https://hub.docker.com/repository/docker/andrewrlapointe/bvbrc/general"
    }

    input {
        File read1
        File read2
        String sample_name
        String user_name
        String password
    }

    runtime {
        docker: 'andrewrlapointe/bvbrc:latest'

    }

    command <<<
        python3 /bin/bvbrc_login.py ~{user_name} ~{password}
        python3 /bin/genome_asm.py ~{user_name} ~{sample_name} ~{read1} ~{read2}
    >>>

    output {
        String testOut = "TODO"
    }
}

task run_genome_annotation {
    meta {
        author: ""
        email: ""
        description: "Task for running the BV-BRC genome annotation tool."
        version: "1.0"
    }

    input {
        String user_name
    }

    runtime {
        docker: 'danylmb/bvbrc:1.040'
    }

    command <<<
        echo "Under construction"
    >>>

    output {
        String anno_out = "TODO"
    }
}

task run_genome_analysis {
    meta {
        author: ""
        email: ""
        description: "Task for running the BV-BRC genome analysis tool."
        version: "1.0"
    }

    input {
        String user_name
    }

    runtime {
        docker: 'danylmb/bvbrc:1.040'
    }

    command <<<
        echo "Under construction"
    >>>

    output {
        String amr_out = "TODO"
    }
}