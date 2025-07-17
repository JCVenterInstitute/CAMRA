version 1.0

import "../Tasks/task_amrfinder.wdl" as amrfinder
import "../Tasks/task_abricate.wdl" as abricate
import "../Tasks/task_hamronization.wdl" as hamronize
import "../Tasks/task_resfinder.wdl" as resfinder
import "../Tasks/task_utilities.wdl" as utilities
import "../Tasks/task_rgi.wdl" as rgi
import "../Tasks/task_bvbrc.wdl" as bvbrc


workflow amr_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Analysis on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        File assembly
         String sample_name
        String organism
 
    }

    # Task to combine genus and species


    call abricate.run_Abricate{
        input:
            assembly = assembly,
            sample_name = sample_name
    }




    output {
    
        # Abricate
        # File abricate_ncbiDB_tsv_output = run_Abricate.abricate_ncbiDB_tsv_output
        File abricate_cardDB_tsv_output = run_Abricate.abricate_cardDB_tsv_output
        # File abricate_resfinderDB_tsv_output = run_Abricate.abricate_resfinderDB_tsv_output
        File abricate_vfdb_tsv_output = run_Abricate.abricate_vfdb_tsv_output
        File abricate_argannotDB_tsv_output = run_Abricate.abricate_argannotDB_tsv_output

        File abricate_DB_version = run_Abricate.abricate_DB_version
        String abricate_version = run_Abricate.abricate_version
        String abricate_date = run_Abricate.abricate_date
    }

}

