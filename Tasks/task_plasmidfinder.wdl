version 1.0

task run_PlasmidFinder {
    meta {
    description: "PlasmidFinder is a tool for the identification and typing of Plasmid Replicons in Whole-Genome Sequencing (WGS)."
    docker:"https://hub.docker.com/layers/staphb/plasmidfinder/2.1.6/images/sha256-caa660124f6fbb0f881b9de6f174b24a6277205e4d734b1cb4b831121b2c769a?context=explore"
    cite: "Carattoli, A., Hasman, H. (2020). PlasmidFinder and In Silico pMLST: Identification and Typing of Plasmid Replicons in Whole-Genome Sequencing (WGS). In: de la Cruz, F. (eds) Horizontal Gene Transfer. Methods in Molecular Biology, vol 2075. Humana, New York, NY. https://doi.org/10.1007/978-1-4939-9877-7_20"
    } 

    runtime{
        docker: 'staphb/plasmidfinder:2.1.6'
    }

    input {
        File assembly
        String sample_name
        File database 
        #these people have a db as a string! 
        #String database # https://github.com/theiagen/public_health_bacterial_genomics/blob/77265020e141ded2f3a5e2af7fe55afd13167abc/tasks/gene_typing/task_plasmidfinder.wdl#L30
    }
    command <<<
        # plasmid finder version is from the docker container
        echo "2.1.6" | tee VERSION 

        # plasmid finder database 
        echo "$(basename ~{database})" | tee DB_VERSION
        date | tee DATE

        #Checking if assembly is zipped
        if [[ "~{assembly}" == *.gz ]]; then
            echo "We need to unzip the assembly."
            # Uncompress the file
            gunzip -c ~{assembly} > assembly.fasta && echo "    > unzip completed"
        else 
            echo "We do not need to unzip."
            mv ~{assembly} assembly.fasta && echo "    > mv successful"
        fi

        #Checking if database is zipped
        if [[ "~{database}" == *.gz ]]; then
            # Uncompress the database
            echo "We need to unzip the database."
            gunzip -c ~{database} > plasmidfinder_db && echo "    > unzip completed"
        else 
            echo "We do not need to unzip."
            mv ~{database} plasmidfinder_db && echo "    > mv successful"
        fi

        #echo "making directory" && mkdir "~{sample_name}_plasmidfinder_results" && echo "directyory made"
        #echo "running plasmid finder" && plasmidfinder.py -i assembly.fasta -p plasmidfinder_db && echo "    > plasmidfiner good"
        mkdir plasmidfinder_output
        echo "running plasmid finder" && plasmidfinder.py -i assembly.fasta -x -o plasmidfinder_output && echo "    > plasmidfiner good" #this works

        #cat data.json | jq ".plasmidfinder.results" | tee PLASMIDFINDER_RESULTS.json
        #grep '"plasmidfinder.results"' data.json | tee PLASMIDFINDER_RESULTS.json
        grep -o -i 'hit_id' plasmidfinder_output/data.json | wc -l | tee HIT_QUANTITY
        tail -n +2 plasmidfinder_output/results_tab.tsv | cut -f2 | tr '\n' ', ' | tee PLASMID_LIST
        rm -rd plasmidfinder_db
        rm assembly.fasta

    >>>

    output {

        File plasmidfinder_tsv_output = "plasmidfinder_output/results_tab.tsv"
        File plasmidfinder_seq_output = "plasmidfinder_output/Hit_in_genome_seq.fsa"
        String plasmidfinder_plasmids_list = read_string("PLASMID_LIST")
        String plasmidfinder_qty_hits = read_string("HIT_QUANTITY")
        String plasmidfinder_date = read_string("DATE")
        String plasmidfinder_version = read_string("VERSION")
        String plasmidfinder_db_version = read_string("DB_VERSION")
    }

}
