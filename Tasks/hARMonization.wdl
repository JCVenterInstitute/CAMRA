version 1.0

task run_Hamronize {
    meta {
    description: "Combine the outputs of disparate antimicrobial resistance gene detection tools into a single unified format."
    gitrepository: "https://github.com/pha4ge/hAMRonization"
    docker:"https://hub.docker.com/r/finlaymaguire/hamronization"
    }

    input {
        Array[File] AMR_files
        Array[File] VIR_files
    }
    runtime{
        docker: 'danylmb/hamronize:v1.1.4-2'
    }

    command <<<
        mkdir AMR_hAMRonization
        for amr_file in ~{sep=" " AMR_files}; do
        
            amr_name=$(basename "$amr_file")
            program="${amr_name%%_*}"
            echo $program
            if [[ $program == "abricate" ]]; then
                hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 1.0.1 --reference_database_version 1.0.0  $amr_file
            fi
            if [[ $program == "amrfinderplus" ]]; then
                hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 3.12.8 --reference_database_version 1.0.0 --input_file_name $amr_file $amr_file
            fi
            if [[ $program == "resfinder" ]]; then 
                hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 4.5.0 --reference_database_version 2.3.1 --input_file_name $amr_file $amr_file
            fi
        done
        hamronize summarize -o hamronize_amr_output.tsv -t tsv AMR_hAMRonization/*

        mkdir VIR_hAMRonization
        for vir_file in ~{sep=" " VIR_files}; do
        
            vir_name=$(basename "$vir_file")
            program="${vir_name%%_*}"
            echo $program
            if [[ $program == "abricate" ]]; then
                hamronize $program --format tsv --output VIR_hAMRonization/"H-$vir_name" --analysis_software_version 1.0.1 --reference_database_version 1.0.0  $vir_file
            fi
            if [[ $program == "amrfinderplus" ]]; then
                hamronize $program --format tsv --output VIR_hAMRonization/"H-$vir_name" --analysis_software_version 3.12.8 --reference_database_version 1.0.0 --input_file_name $vir_file $vir_file
            fi
            if [[ $program == "resfinder" ]]; then 
                hamronize $program --format tsv --output VIR_hAMRonization/"H-$vir_name" --analysis_software_version 4.5.0 --reference_database_version 2.3.1 --input_file_name $vir_file $vir_file
            fi
        done
        hamronize summarize -o hamronize_vir_output.tsv -t tsv VIR_hAMRonization/*



        >>>

    output{
        File hAMRonization_amr_output = "hamronize_amr_output.tsv"
        File hAMRonization_vir_output = "hamronize_vir_output.tsv"
    }

}