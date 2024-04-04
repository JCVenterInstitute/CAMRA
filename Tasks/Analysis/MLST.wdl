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
        String sample_name
        String organism
        String? docker = 'staphb/mlst:2.23.0-2024-01'
        String? scheme
        Float? minid
        Float? mincov
        Float? minscore

    }
    runtime{
        docker: docker
    }

    command <<<

    echo $(mlst --version 2>&1) | sed 's/mlst //' | tee VERSION
    # Update the database
        #Q? would this be useful? it seems to take some time to create this database.
        # the MLST tool author updates his db every so often, but it is quite different when i update it
    /mlst-2.23.0/scripts/mlst-download_pub_mlst -d pubmlst | bash
    # renames the old database
    mv $(find ../mlst*/db -name pubmlst) "$(find ../mlst*/db -name pubmlst).old"
    # moves the new database
    mv pubmlst/ ../mlst*/db/


    #Q? i have multiple places where i could get taxa from lab, mash, kraken, mlst, how do i know which one to use? 
    mlst --novel ~{sample_name}_novel_mlst.fasta --csv ~{assembly} > ~{sample_name}-MLST-NoScheme.csv 
    mlst \
    ~{'--scheme ' + scheme} \
    ~{'--minid ' + minid} \
    ~{'--mincov ' + mincov} \
    ~{'--minscore ' + minscore} \
    --novel ~{sample_name}_novel_mlst_alleles.fasta \
    ~{assembly} \
    >> ~{sample_name}_ts_mlst.tsv


    >>>
    output{
    }
}
