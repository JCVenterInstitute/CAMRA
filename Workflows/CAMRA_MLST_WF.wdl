version 1.0

import "../Tasks/task_mlst.wdl" as mlst

workflow get_mlst {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run QC on a SINGLE pared-end illumina Assembly. Runs Mash, CheckM (inputs mash taxa), Merqury."
    }

    parameter_meta {
        sample_name     :   "Name of sample isolate"
        assembly        :   "Assembly in fasta or fasta.gz format"
        pubmlst_DB      :   "Optional pubmlst_DB for mlst mapping"
    }

    input {
        String sample_name
        File assembly
        File? pubmlst_DB    
    }
    
    call mlst.run_MLST if !defined(pubmlst_DB) {
        input:
            assembly = assembly
    }

    call mlst.run_MLST_pubmlst_DB if defined(pubmlst_DB) {
        input:
            assembly = assembly,
            pubmlst_DB = pubmlst_DB
    }

    output {

        #MLST
        String? tsMLST_version = mlst.run_MLST.tsMLST_version
        String? tsMLST_date = mlst.run_MLST.tsMLST_date
        File? tsMLST_tsv_output = mlst.run_MLST.tsMLST_tsv_output
        String? tsMLST_scheme = mlst.run_MLST.tsMLST_scheme 
        String? tsMLST_seqtype = mlst.run_MLST.tsMLST_seqtype 
        String? tsMLST_alleles = mlst.run_MLST.tsMLST_alleles

        #MLST + pubmlst_DB
        String? tsMLSTdb_version = mlst.run_MLST_pubmlst_DB.tsMLST_version
        String? tsMLSTdb_date = mlst.run_MLST_pubmlst_DB.tsMLST_date
        File? tsMLSTdb_tsv_output = mlst.run_MLST_pubmlst_DB.tsMLST_tsv_output
        String? tsMLSTdb_scheme = mlst.run_MLST_pubmlst_DB.tsMLST_scheme 
        String? tsMLSTdb_seqtype = mlst.run_MLST_pubmlst_DB.tsMLST_seqtype 
        String? tsMLSTdb_alleles = mlst.run_MLST_pubmlst_DB.tsMLST_alleles
    }


}
