version 1.0


task run_RGI { 
    meta{
        description: "This application is used to predict antibiotic resistome(s) from protein or nucleotide data based on homology and SNP models. "
        repository:"https://github.com/arpcard/rgi.git"
        cite:"Alcock et al. 2023. CARD 2023: Expanded Curation, Support for Machine Learning, and Resistome Prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, D690-D699."
    }

    input {
        File assembly
        Boolean ? diamond = false
        Boolean ? blast = true
    }
    runtime {
        docker:"quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0"
        cpu: 4
        memory: "16GB"
        disks:  "local-disk " + 50 + " HDD"
        disk: 50 + " GB" # TES
    }


    command <<<
        rgi main --version | tee RGI_VERSION
        rgi database --version | tee CARD_DB_VERSION
        date | tee DATE

        if [ ~{blast} == true ]; then
            rgi main \
                -i ~{assembly} \
                -o rgi_blast_report \
                -t contig \
                -a BLAST \
                -d wgs \
                --split_prodigal_jobs \
                --clean 
        fi

        if [[ ~{diamond} == true || ~{blast} == false ]]; then
            rgi main \
                -i ~{assembly} \
                -o rgi_diamond_report \
                -t contig \
                -a DIAMOND \
                -d wgs \
                --split_prodigal_jobs \
                --clean
        fi

    >>>
    output {
        String rgi_CARD_DB_version = read_string("CARD_DB_VERSION")
        String rgi_version = read_string("RGI_VERSION")
        String rgi_date = read_string("DATE")

        File? rgi_CARD_diamond_tsv_output = "rgi_diamond_report.txt"
        File? rgi_CARD_blast_tsv_output = "rgi_blast_report.txt"
        File? rgi_CARD_diamond_json_output = "rgi_diamond_report.json"
        File? rgi_CARD_blast_json_output = "rgi_blast_report.json"
    }
    

}