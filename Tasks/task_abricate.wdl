version 1.0

task run_Abricate {
    meta {
    description: "Mass screening of contigs for antimicrobial resistance or virulence genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB."
    gitrepository: "https://github.com/tseemann/abricate"
    docker:"https://hub.docker.com/r/staphb/abricate"
    cite: "Seemann T, Abricate, Github https://github.com/tseemann/abricate; This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust"
    } 

    input {
        File assembly
        String sample_name
        #File pubmlst_DB
    }
    runtime{
        docker: 'staphb/abricate:1.0.1-insaflu-220727'
    }

    command <<<

        date | tee DATE
        echo $(abricate --version 2>&1) | sed 's/abricate //' | tee VERSION

        #Updating the databases
        abricate-get_db --db ncbi --force
        abricate-get_db --db card --force
        abricate-get_db --db resfinder --force
        #abricate-get_db --db megares --force #can not update
        abricate-get_db --db vfdb --force
        abricate-get_db --db argannot --force

        echo $(abricate --list) | tee DB_VERSION
        grep -w -f <(echo -e "ncbi\ncard\nresfinder\nvfdb\nargannot") DB_VERSION > DB_VERSION

        # Function to run abricate for a specific database and save the hits
        run_abricate() {
            local db="$1"
            local output_file="abricate_~{sample_name}_${db}_hits.tsv"
            abricate --threads 4 --db "$db" ~{assembly} > "$output_file"

            # Parse out gene names into list of strings, comma-separated, final comma at end removed by sed
            local genes=$(awk -F '\t' '{ print $6 }' "$output_file" | tail -n+2 | tr '\n' ',' | sed 's/.$//')

            # If variable for list of genes is EMPTY, write string saying it is empty
            if [ -z "$genes" ]; then
                genes="No genes detected by ABRicate $db db"
            fi

            # Output genes to file
            echo "$genes" > "ABRICATE_GENES_${db^^}"
        }

        # Run abricate for each database
        databases=("ncbi" "card" "resfinder" "vfdb" "argannot")
        for db in "${databases[@]}"; do
            run_abricate "$db"
        done
        >>>

    output{
        File abricate_DB_version = "DB_VERSION"
        String abricate_version = read_string("VERSION")
        String abricate_date = read_string("DATE")

        File abricate_ncbiDB_tsv_output = "abricate_~{sample_name}_ncbi_hits.tsv"
        File abricate_cardDB_tsv_output = "abricate_~{sample_name}_card_hits.tsv"
        File abricate_resfinderDB_tsv_output = "abricate_~{sample_name}_resfinder_hits.tsv"
        File abricate_vfdb_tsv_output = "abricate_~{sample_name}_vfdb_hits.tsv"
        File abricate_argannotDB_tsv_output = "abricate_~{sample_name}_argannot_hits.tsv"

    }





}