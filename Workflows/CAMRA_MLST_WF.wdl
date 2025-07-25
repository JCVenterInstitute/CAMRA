version 1.0

import "../Tasks/task_mlst.wdl" as mlst

Workflow assembly_qc {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run QC on a SINGLE pared-end illumina Assembly. Runs Mash, CheckM (inputs mash taxa), Merqury."
    }

    parameter_meta {
        sample_name     :   "Name of sample isolate"
        read1           :   "Raw read 1 fastq.gz or fastq file"
        read2           :   "Raw read 2 fastq.gz or fastq file"
        assembly        :   "Assembly in fasta or fasta.gz format"
        pubmlst_DB      :   "Optional pubmlst_DB for mlst mapping"
    }

    input{
        String sample_name
        File read1
        File read2
        File assembly
        File pubmlst_DB = "None"
    }
    
    if ( "~{pubmlst_DB}" == "None") {
        call mlst.run_MLST {
            input:
                assembly = assembly
        }
    }
    
    if ( "~{pubmlst_DB}" != "None" ) {
        call mlst.run_MLST_pubmlst_DB {
            input: 
                assembly = assembly,
                pubmlst_DB = pubmlst_DB
        }
    } 


    output{

        #MLST
        String? tsMLST_version = run_MLST.tsMLST_version
        String? tsMLST_date = run_MLST.tsMLST_date
        File? tsMLST_tsv_output = run_MLST.tsMLST_tsv_output
        String? tsMLST_scheme = run_MLST.tsMLST_scheme 
        String? tsMLST_seqtype = run_MLST.tsMLST_seqtype 
        String? tsMLST_alleles = run_MLST.tsMLST_alleles

        #MLST + pubmlst_DB
        String? tsMLSTdb_version = run_MLST_pubmlst_DB.tsMLST_version
        String? tsMLSTdb_date = run_MLST_pubmlst_DB.tsMLST_date
        File? tsMLSTdb_tsv_output = run_MLST_pubmlst_DB.tsMLST_tsv_output
        String? tsMLSTdb_scheme = run_MLST_pubmlst_DB.tsMLST_scheme 
        String? tsMLSTdb_seqtype = run_MLST_pubmlst_DB.tsMLST_seqtype 
        String? tsMLSTdb_alleles = run_MLST_pubmlst_DB.tsMLST_alleles
    }


}
