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
        docker: 'danylmb/hamronize:v1.1.4'
    }

    command <<<
        mkdir hAMRonization

        for amr_file in ~{sep=" " AMR_files}; do
            amr_name=$(basename "$amr_file")
            program="${amr_name%%_*}"
            if [[ $program == "abricate" ]]; then
                hamronize $program --format tsv --output hAMRonization/"H-$amr_name" --analysis_software_version 1.0.1 --reference_database_version 1.0.0  $amr_file
            fi
            if [[ $program == "amrfinderplus" ]]; then
                hamronize $program --format tsv --output hAMRonization/"H-$amr_name" --analysis_software_version 3.12.8 --reference_database_version 1.0.0 --input_file_name $amr_file $amr_file
            fi

            hamronize summarize -o hamronize_output.tsv -t tsv hAMRonization/*
        done


        





        >>>

    output{

    }

}