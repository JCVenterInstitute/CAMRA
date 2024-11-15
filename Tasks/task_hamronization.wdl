version 1.0

task run_hAMRonize {
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
        docker: 'danylmb/hamronize:v1.1.4-build10'
    }

    command <<<

        date | tee DATE
        hamronize --version | tee VERSION

        ##################################################
        # FUNCTIONS
        ##################################################
        check_dataframe_rows() {
            local file_path="$1"  # Take the first argument as the input (file path)

            # Check if the file exists and is not empty
            if [[ ! -s "$file_path" ]]; then
                echo "    No data available in the DataFrame or file does not exist."
                return 1  # Return with an error code to indicate no data or file not found
            fi

            # Count the number of lines excluding the header
            local row_count=$(tail -n +2 "$file_path" | wc -l)

            # Check if there are rows present
            if [[ $row_count -gt 0 ]]; then
                echo "    The DataFrame has rows."
                return 0  # Return success code
            else
                echo "    The DataFrame is empty (only headers)."
                return 1  # Return with an error code to indicate no rows
            fi
        }

        echo "##################################################"
        echo "START AMR HARMONIZATION "
        echo "##################################################"
        mkdir AMR_hAMRonization

        # Iterate thought all available AMR files and standarize them.
        for amr_file in ~{sep=" " AMR_files}; do
            amr_name=$(basename "$amr_file") # Get the File's basename
            program="${amr_name%%_*}" # Get the tool the used to create the file (eg abricate, amrfinder, resfinder)
            echo $program

            if [[ $program == "abricate" ]]; then
                echo "    $amr_file = abricate"
                if check_dataframe_rows $amr_file; then
                    echo "    Starting hamronization of $amr_name"
                    hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 1.0.1 --reference_database_version 1.0.0  $amr_file
                fi   
            fi

            if [[ $program == "amrfinderplus" ]]; then
                echo "    $amr_file = amrfinderplus"
                if check_dataframe_rows $amr_file; then
                    echo "    Starting hamronization of $amr_name"
                    hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 3.12.8 --reference_database_version 1.0.0 --input_file_name $amr_file $amr_file
                fi
            fi
            if [[ $program == "resfinder" ]]; then 
                echo "    $amr_file = resfinder"
                if check_dataframe_rows $amr_file; then
                    echo "    Starting hamronization of $amr_name"
                    hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 4.5.0 --reference_database_version 2.3.1 --input_file_name $amr_file $amr_file
                fi
            fi

            if [[ $program == "rgi" ]]; then 
                echo "    $amr_file = rgi"
                if check_dataframe_rows $amr_file; then
                    echo "    Starting hamronization of $amr_name"
                    hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --analysis_software_version 6.0.3 --reference_database_version 3.3.0 --input_file_name $amr_file $amr_file
                fi
            fi

            if [[ $program == "bvbrc" ]]; then 
                echo "    $amr_file = bvbrc"
                if check_dataframe_rows $amr_file; then
                    echo "    Starting hamronization of $amr_name"
                    hamronize $program --format tsv --output AMR_hAMRonization/"H-$amr_name" --input_file_name $amr_file $amr_file
                fi
            fi
        done

        # Check if there are any files in the directory
        if [ -n "$(ls -A AMR_hAMRonization/ 2>/dev/null)" ]; then
            # If the directory is not empty, run the hamronize command
            hamronize summarize -o hamronize_amr_output.tsv -t tsv AMR_hAMRonization/* && echo "    AMR Harmonisation DONE"
        else
            # If the directory is empty, print a message
            echo "The directory AMR_hAMRonization/ is empty. "
            touch hamronize_amr_output.tsv
        fi

        echo "##################################################"
        echo "START VIRULENCE HARMONIZATION "
        echo "##################################################"
        mkdir VIR_hAMRonization

        for vir_file in ~{sep=" " VIR_files}; do
            vir_name=$(basename "$vir_file")
            program="${vir_name%%_*}"
            echo $program

            if [[ $program == "abricate" ]]; then
                if check_dataframe_rows "$program"; then
                    echo "    Starting hamronization of $vir_name"
                    hamronize $program --format tsv --output VIR_hAMRonization/"H-$vir_name" --analysis_software_version 1.0.1 --reference_database_version 1.0.0  $vir_file
                fi
            fi

            if [[ $program == "amrfinderplus" ]]; then
                if check_dataframe_rows "$program"; then
                    hamronize $program --format tsv --output VIR_hAMRonization/"H-$vir_name" --analysis_software_version 3.12.8 --reference_database_version 1.0.0 --input_file_name $vir_file $vir_file
                fi
            fi

            if [[ $program == "resfinder" ]]; then 
                if check_dataframe_rows "$program"; then
                    hamronize $program --format tsv --output VIR_hAMRonization/"H-$vir_name" --analysis_software_version 4.5.0 --reference_database_version 2.3.1 --input_file_name $vir_file $vir_file
                fi
            fi
        done

        # Check if there are any files in the directory
        if [ -n "$(ls -A VIR_hAMRonization/ 2>/dev/null)" ]; then
            # If the directory is not empty, run the hamronize command
            hamronize summarize -o hamronize_vir_output.tsv -t tsv VIR_hAMRonization/* && echo "    Virulence Harmonisation DONE"
        else
            # If the directory is empty, print a message
            echo "The directory AMR_hAMRonization/ is empty. "
            touch hamronize_vir_output.tsv
        fi


        echo "##################################################"
        echo "START AMR TERM CONSOLIDATION "
        echo "##################################################"


        if check_dataframe_rows hamronize_amr_output.tsv; then # Check 
            python3 /usr/bin/amr-term-consolidation.py hamronize_amr_output.tsv && echo "    Consolidation DONE"
        else
            echo "    No hamronize_amr_output.tsv"
            touch consolidation_isna.tsv consolidation_all.tsv consolidation_amr_over98identity.tsv consolidation_amr_allidentity.tsv
        fi

        >>>

    output{
        String hAMRonization_version = read_string("VERSION")
        String hAMRonization_date = read_string("DATE")

        
        File hAMRonization_amr_output = "hamronize_amr_output.tsv"
        File hAMRonization_vir_output = "hamronize_vir_output.tsv"
        String amrtermconsolidation_version = "1.0.0"
        File amrtermconsolidation_isna = "consolidation_isna.tsv"
        File amrtermconsolidation_all = "consolidation_all.tsv"
        File amrtermconsolidation_over98 = "consolidation_amr_over98identity.tsv"
        File amrtermconsolidation_allidentity = "consolidation_amr_allidentity.tsv"
    }

}