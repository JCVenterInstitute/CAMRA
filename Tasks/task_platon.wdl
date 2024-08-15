version 1.0

task run_platon {
    meta {
        description: ""
        version: "1.0"
        dockerhub: ""
    }

    runtime {
        docker: "andrewrlapointe/platon:1.0"
    }

    input {
        File assembly
    }

    command <<<
        echo $PATH
        platon --db /db ~{assembly}
    >>>

    output {
        String mob_out = "TODO"
    }
}