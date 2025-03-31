version 1.0

task run_hAMRonize {
    meta {
    description: "Combine the outputs of disparate antimicrobial resistance gene detection tools into a single unified format."
    gitrepository: "https://github.com/pha4ge/hAMRonization"
    docker:"https://hub.docker.com/r/finlaymaguire/hamronization"
    }

    input {

        #AMR Output 
        # File abricate_ncbiDB_tsv_output
        File abricate_cardDB_tsv_output
        # File abricate_resfinderDB_tsv_output
        File abricate_argannotDB_tsv_output 

        File amrfinder_amr_output

        File resfider_asm_output
        File? resfinder_read_output

        File? rgi_CARD_diamond_tsv_output
        File? rgi_CARD_blast_tsv_output 

        File? bvbrc_amr_file

        # Virulence Output

        Array[File]+ VIR_files

        # CARD DB
        File CARD_aro_obo
        File CARD_protein_homolog
        File CARD_protein_variant
        File CARD_nucleotide_homolog
        File CARD_nucleotide_variant

        # Assembly
        File assembly
    }
    runtime{
        docker: 'danylmb/hamronize:v1.1.4-build17'
        continueOnReturnCode: [0, 1]
    }

    command <<<

        date | tee DATE
        hamronize --version | tee VERSION

        echo "/WHERE AM I"
        echo "$(pwd)"

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

        # Initialize an empty array
        file_paths=()

        # Mandatory files
        file_paths+=( ~{abricate_cardDB_tsv_output} ~{abricate_argannotDB_tsv_output} ~{amrfinder_amr_output} ~{resfider_asm_output})

        # Conditionally add optional files if they are defined
        if ~{defined(bvbrc_amr_file)}; then
            file_paths+=(~{bvbrc_amr_file})
        fi

        if ~{defined(resfinder_read_output)}; then
            file_paths+=(~{resfinder_read_output})
        fi

        if ~{defined(rgi_CARD_blast_tsv_output)}; then
            file_paths+=(~{rgi_CARD_blast_tsv_output})
        fi

        if ~{defined(rgi_CARD_diamond_tsv_output)}; then
            file_paths+=(~{rgi_CARD_diamond_tsv_output})
        fi


        # Iterate thought all available AMR files and standarize them.
        for amr_file in ${file_paths[@]}; do
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
            echo "The directory AMR_hAMRonization/ is not empty. "
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
        # cd /data
        mkdir AMR_Term_Consolidation

        # Check if file exists and is not empty
        if [[ -f "hamronize_amr_output.tsv" && -s "hamronize_amr_output.tsv" ]]; then
            if [[ "~{assembly}" == *.fasta.gz || "~{assembly}" == *.fa.gz ]]; then 
                gunzip -c ~{assembly} > assembly.fasta 
            elif [[ "~{assembly}" == *.fasta || "~{assembly}" == *.fa ]]; then
                mv ~{assembly} assembly.fasta
            fi
        
            python3 /opt/AMR_Term_Consolidation/AMR_Term_Consolidation.py \
                "./hamronize_amr_output.tsv" \
                ~{CARD_aro_obo} \
                assembly.fasta \
                ~{CARD_protein_homolog} \
                ~{CARD_protein_variant} \
                ~{CARD_nucleotide_homolog} \
                ~{CARD_nucleotide_variant}
        else
            echo "    No hamronize_amr_output.tsv"
            touch consolidation_isna.tsv consolidation_all.tsv consolidation_amr_over98identity.tsv consolidation_amr_allidentity.tsv
        fi

        rm *.ntf *.nto *.ndb *.nhr *.nin *.njs *.not *.nsq *.pdb *.phr *.pin *.pjs *.pot *.psq *.ptf *.pto assembly.fasta

    >>>

    output{
        String hAMRonization_version = read_string("VERSION")
        String hAMRonization_date = read_string("DATE")

        
        File? hAMRonization_amr_output = "hamronize_amr_output.tsv"
        File? hAMRonization_vir_output = "hamronize_vir_output.tsv"

        File amrtermconsolidation_log = "term_consolidation.log"

        File hAMRonization_HARMONIZED_TERMS = "HARMONIZED_TERMS.tsv"
        File hAMRonization_CONSOLIDATED_TERMS = "CONSOLIDATED_TERMS.tsv"


        String amrtermconsolidation_version = "1.0.0"

        

        File? amrtermconsolidation_isna = "consolidation_isna.tsv"
        File? amrtermconsolidation_all = "consolidation_all.tsv"
        File? amrtermconsolidation_over98 = "consolidation_amr_over98identity.tsv"
        File? amrtermconsolidation_allidentity = "consolidation_amr_allidentity.tsv"
    }

}