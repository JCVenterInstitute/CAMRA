version 1.0

task run_Hamronize {
    meta {
    description: "Combine the outputs of disparate antimicrobial resistance gene detection tools into a single unified format."
    gitrepository: "https://github.com/pha4ge/hAMRonization"
    docker:"https://hub.docker.com/r/finlaymaguire/hamronization"
    }

    input {
        Array[File] AMR_files
        String sample_name
    }
    runtime{
        docker: 'finlaymaguire/hamronization:1.0.0'
    }

    command <<<
        hamronize --help
        echo ~{sep=", "AMR_files} && echo "it worked"

        >>>

    output{

    }

}