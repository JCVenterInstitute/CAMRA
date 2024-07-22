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
        docker: 'andrewrlapointe/bvbrc:3.0'

    }

    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py assembly ~{username} ~{sample_name} ~{read1} ~{read2}
    >>>

    output {
        String testOut = "TODO"
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
        docker: 'andrewrlapointe/bvbrc:3.0'
    }

    # output path could be changed to be relative to the contigs file location to reduce the number of inputs
    command <<<
        python3 /bin/bvbrc_login.py ~{username} ~{password}
        python3 /bin/bvbrc_jobs.py cga ~{contigs_file} ~{output_path} ~{sample_name}_output ~{scientific_name} ~{taxonomy_id}
    >>>

    output {
        String amr_out = "TODO"
    }
}