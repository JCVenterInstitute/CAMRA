version 1.0

import "../Tasks/task_rgi.wdl" as rgi

workflow amr_analysis   {
    input {
        File assembly
        Boolean update_CARD = false
    }
    call rgi.run_RGI {
        input:
            assembly = assembly,
            update_CARD = update_CARD 
    }

    output{
    # RGI
        String rgi_CARD_DB_version = run_RGI.rgi_CARD_DB_version
        String rgi_version = run_RGI.rgi_version
        String rgi_date = run_RGI.rgi_date

        File rgi_CARD_diamond_tsv_output = run_RGI.rgi_CARD_diamond_tsv_output
        File rgi_CARD_blast_tsv_output = run_RGI.rgi_CARD_blast_tsv_output
        File rgi_CARD_diamond_json_output = run_RGI.rgi_CARD_diamond_json_output
        File rgi_CARD_blast_json_output = run_RGI.rgi_CARD_blast_json_output
    }  
}
