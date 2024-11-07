version 1.0

task run_mob_suite {
    meta {
        description: ""
        version: "1.0"
        dockerhub: ""
    }

    input {
        File assembly
    }

    runtime {
        docker: "kbessonov/mob_suite:3.0.3"
    }

    command <<<
        echo TODO
    >>>

    output {
        String mob_out = "TODO"
    }
}