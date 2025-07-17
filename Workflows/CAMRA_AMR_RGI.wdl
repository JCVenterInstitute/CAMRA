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


    call rgi.run_RGI {
        input:
            assembly = assembly
    }

   


    output {
        # RGI
        String rgi_CARD_DB_version = run_RGI.rgi_CARD_DB_version
        String rgi_version = run_RGI.rgi_version
        String rgi_date = run_RGI.rgi_date

        File? rgi_CARD_diamond_tsv_output = run_RGI.rgi_CARD_diamond_tsv_output
        File? rgi_CARD_blast_tsv_output = run_RGI.rgi_CARD_blast_tsv_output
        File? rgi_CARD_diamond_json_output = run_RGI.rgi_CARD_diamond_json_output
        File? rgi_CARD_blast_json_output = run_RGI.rgi_CARD_blast_json_output
    }

}

