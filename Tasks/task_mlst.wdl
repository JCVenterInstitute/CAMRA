version 1.0

task run_MLST {
    meta {
    description: "Torsten Seeman's (TS) automatic MLST scans contig files against traditional PubMLST typing schemes"
    gitrepository: "https://github.com/tseemann/mlst"
    docker:"https://hub.docker.com/r/staphb/mlst"
    cite: "Seemann T, mlst Github https://github.com/tseemann/mlst; This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust"
    } 


    input {
        File assembly
    }
    runtime{
        docker: 'staphb/mlst:2.23.0-2024-01'
    }

    command <<<

        echo $(mlst --version 2>&1) | sed 's/mlst //' | tee VERSION
        date | tee DATE
        echo "this is the assembly info $(ls -l ~{assembly})"

        echo "Running MLST" && \
        mlst \
        --novel \
        novel_mlst.fasta \
        ~{assembly} \
        >> ts_mlst_output.tsv  && cd \
        echo "MLST done"

        cut -f2 "ts_mlst_output.tsv" | tail -n 1 | tee MLST_SCHEME
        cut -f3 "ts_mlst_output.tsv" | tail -n 1 | tee MLST_SEQTYPE
        cut -f 4- "ts_mlst_output.tsv" | tail -n 1 | sed -e 's|\t|,|g' | tee MLST_ALLELICPROFILE

    >>>

    output{
        String tsMLST_version = read_string("VERSION")
        String tsMLST_date = read_string("DATE")

        File tsMLST_tsv_output = "ts_mlst_output.tsv"

        String tsMLST_scheme = read_string("MLST_SCHEME")
        String tsMLST_seqtype = read_string("MLST_SEQTYPE")
        String tsMLST_alleles = read_string("MLST_ALLELICPROFILE")

    }
}

task run_MLST_pubmlst_DB {
    meta {
    description: "Torsten Seeman's (TS) automatic MLST scans contig files against traditional PubMLST typing schemes"
    gitrepository: "https://github.com/tseemann/mlst"
    docker:"https://hub.docker.com/r/staphb/mlst"
    cite: "Seemann T, mlst Github https://github.com/tseemann/mlst; This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust"
    } 


    input {
        File assembly
        File pubmlst_DB
    }
    
    runtime{
        docker: 'staphb/mlst:2.23.0-2024-01'
    }

    command <<<

        echo $(mlst --version 2>&1) | sed 's/mlst //' | tee VERSION
        date | tee DATE
        echo "this is the assembly info $(ls -l ~{assembly})"

        echo "Pubmlst was input."

        # Remove the pubmlst database belonging to the docker
        rm -rd $(find /mlst*/db -name pubmlst)

        # Make a new pubmlst directory, open the input pubmlst_DB tar.gz file into the new pubmlst directory
        mkdir "$(find / -type d -name mlst*)/db/pubmlst" && echo "new pubmlst made"
        tar -xzf ~{pubmlst_DB} -C "$(find / -type d -name mlst*)/db/pubmlst" --strip-components=1 && echo "input pubmlst was unziped"

        # Records the location where the wdl is executing
        original_pwd=$(pwd)

        # Enter the directory where the MLST scripts are located
        # Run mlst-make_blast_db script. This will generatre the new blast_db from the input pubmlst_DB
        cd $(dirname $(find / -type f -name "mlst-make_blast_db")) && ./mlst-make_blast_db && echo "mlst-make_blast_db was ran"

        # Return to the wdl rexecution directory
        cd $original_pwd

        #MLST will run
        echo "Running MLST"
        mlst \
        --novel \
        novel_mlst.fasta \
        ~{assembly} \
        >> "ts_mlst_output.tsv"
        
        #Parse the output of st_mlst.tsv file into their own files.
        cut -f2 "ts_mlst_output.tsv" | tail -n 1 | tee MLST_SCHEME
        cut -f3 "ts_mlst_output.tsv" | tail -n 1 | tee MLST_SEQTYPE
        cut -f 4- "ts_mlst_output.tsv" | tail -n 1 | sed -e 's|\t|,|g' | tee MLST_ALLELICPROFILE

        >>>

    output{
        String tsMLST_version = read_string("VERSION")
        String tsMLST_date = read_string("DATE")

        File tsMLST_tsv_output = "ts_mlst_output.tsv"

        String tsMLST_scheme = read_string("MLST_SCHEME")
        String tsMLST_seqtype = read_string("MLST_SEQTYPE")
        String tsMLST_alleles = read_string("MLST_ALLELICPROFILE")

    }
}
