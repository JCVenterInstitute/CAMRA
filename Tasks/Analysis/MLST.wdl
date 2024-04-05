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
        #TODO how big is the database?
        # freeze the version of the db in a new task/wf  generate tar file and use as input
        # volatile tag is also an ooption , lets wdl executer know, it will always run, even downstream 
            #See; https://github.com/broadinstitute/viral-pipelines/blob/967ccb9709c5395b074ec3784b1b6af2001c18ba/pipes/WDL/tasks/tasks_utils.wdl#L692
        #broad: they seperate the db and commmands, input are a tar file
    /mlst-2.23.0/scripts/mlst-download_pub_mlst -d pubmlst | bash
    # renames the old database
    mv $(find ../mlst*/db -name pubmlst) "$(find ../mlst*/db -name pubmlst).old"
    # moves the new database
    mv pubmlst/ ../mlst*/db/


    #Q? i have multiple places where i could get taxa from lab, mash, kraken, mlst, how do i know which one to use? 
    # perhaps its good to offer more than one answer, 
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

task get_MLST_db{
    meta {
    description: "Generates the most uptodate pubmlst database"
    gitrepository: "https://github.com/tseemann/mlst"
    docker:"https://hub.docker.com/r/staphb/mlst"
    cite: "Seemann T, mlst Github https://github.com/tseemann/mlst; This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley (Jolley & Maiden 2010, BMC Bioinformatics, 11:595) and sited at the University of Oxford. The development of that website was funded by the Wellcome Trust"
    } 


    input {
        String docker = 'staphb/mlst:2.23.0-2024-01'
    }
    runtime{
        docker: docker
    }
    command <<<
        # Update the database. expect db of ~100M
        /mlst-2.23.0/scripts/mlst-download_pub_mlst -d pubmlst | bash
        db_name="pubmlst-$(date "+%d%m%Y").tar.gz"
        tar -czvf "$db_name" -C pubmlst . 
    >>>
    output{
        File pubmlst_db = "$db_name"
    }
}
