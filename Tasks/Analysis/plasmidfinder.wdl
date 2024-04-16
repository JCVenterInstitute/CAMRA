version 1.0

task run_PlasmidFinder {
    meta {
    description: "PlasmidFinder is a tool for the identification and typing of Plasmid Replicons in Whole-Genome Sequencing (WGS)."
    docker:"https://hub.docker.com/layers/staphb/plasmidfinder/2.1.6/images/sha256-caa660124f6fbb0f881b9de6f174b24a6277205e4d734b1cb4b831121b2c769a?context=explore"
    cite: "Carattoli, A., Hasman, H. (2020). PlasmidFinder and In Silico pMLST: Identification and Typing of Plasmid Replicons in Whole-Genome Sequencing (WGS). In: de la Cruz, F. (eds) Horizontal Gene Transfer. Methods in Molecular Biology, vol 2075. Humana, New York, NY. https://doi.org/10.1007/978-1-4939-9877-7_20"
    } 

    runtime{
        docker: 'staphb/plasmidfinder:2.1.6_2024-03-07'
    }

    input {
        File assembly
        String sample_name
        File database #could be optional. 
        #File pubmlst_DB
    }
    command <<<
        # plasmid finder version is from the docker container
        echo "2.1.6_2024-03-07" | tee VERSION 

        # plasmid finder database 
        echo "$(basename ~{database})" | tee DB_VERSION
        date | tee DATE

        #Checking if assembly is zipped
        if [[ "~{assembly}" == *.gz ]]; then
            # Uncompress the file
            gunzip -c "~{assembly}" > "assembly.fasta"
        else 
            mv ~{assembly} "assembly.fasta"
        if

        #Checking if database is zipped
        if [[ "~{database}" == *.gz ]]; then
            # Uncompress the database
            gunzip -c "~{database}" > "plasmidfinder_db"
        else 
            mv ~{assembly} "plasmidfinder_db"
        if

        mkdir plasmidfinder_results
        plasmidfinder.py -i assembly.fasta -o plasmidfinder_results -p plasmidfinder_db

        cat plasmidfinder_results/data.json | jq ".plasmidfinder.results" | tee PLASMIDFINDER_RESULTS.json
        grep -o -i 'hit_id' plasmidfinder_results/data.json | wc -l | tee HIT_QUANTITY


    >>>

    output {
        File plasmidfinder_json_output = "plasmidfinder_results/data.json"
        File plasmidfinder_results = "PLASMIDFINDER_RESULTS.json"
        String plasmidfinder_qtyhits = read_string("HIT_QUANTITY")
        String plasmidfinder_date = read_string("DATE")
        String plasmidfinder_version = read_string("VERSION")
        String plasmidfinder_dbversion = read_string("DB_VERSION")
    }

}
