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
        docker: 'andrewrlapointe/bvbrc:latest'

    }

    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py assembly ~{username} ~{sample_name} ~{read1} ~{read2}
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
        File contigs_file
        String username
        String password
        String sample_name
        String? scientific_name
    }

    runtime {
        docker: 'andrewrlapointe/bvbrc:latest'
    }

    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py annotation ~{username} ~{sample_name} ~{contigs_file}

    >>>

    output {
        String anno_out = "TODO"
    }
}

task run_genome_analysis {
    meta {
        author: ""
        email: ""
        description: "Task for running the BV-BRC complete genome analysis tool."
        version: "1.0"
    }

    input {
        File contigs_file
        String username
        String password
        String assembly_filepath
        String sample_name
        String? scientific_name
    }

    runtime {
        docker: 'andrewrlapointe/bvbrc:latest'
    }

    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py cga ~{username} ~{sample_name} ~{assembly_filepath}
    >>>

    output {
        String amr_out = "TODO"
    }
}