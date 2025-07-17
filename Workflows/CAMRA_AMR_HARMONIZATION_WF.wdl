version 1.0

import "../Tasks/task_hamronization.wdl" as hamronize


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
        # File? read1
        # File? read2
        File? blast_query
        String sample_name
        String organism
        File amrfinder_amr_output
        File resfider_asm_output
        File? rgi_CARD_diamond_tsv_output
        File rgi_CARD_blast_tsv_output
        File abricate_cardDB_tsv_output
        File abricate_argannotDB_tsv_output
        File abricate_vfdb_tsv_output
        File amrfinder_virulence_output
   }

    # Task to combine genus and species



    call hamronize.run_hAMRonize {
        input:
            assembly = assembly,
            abricate_cardDB_tsv_output  = abricate_cardDB_tsv_output,
            abricate_argannotDB_tsv_output  = abricate_argannotDB_tsv_output,

            amrfinder_amr_output = amrfinder_amr_output,

            resfider_asm_output = resfider_asm_output,
            rgi_CARD_blast_tsv_output = rgi_CARD_blast_tsv_output,

            # Virulence Output

            VIR_files = [
                abricate_vfdb_tsv_output,
                amrfinder_virulence_output
            ]
        }



    output {
            # hAMRonization
        String hAMRonization_version = run_hAMRonize.hAMRonization_version
        String hAMRonization_date = run_hAMRonize.hAMRonization_date
        File? hAMRonization_amr_output = run_hAMRonize.hAMRonization_amr_output
        File? hAMRonization_vir_output = run_hAMRonize.hAMRonization_vir_output
        File hAMRonization_HARMONIZED_TERMS = run_hAMRonize.hAMRonization_HARMONIZED_TERMS
        File hAMRonization_CONSOLIDATED_TERMS = run_hAMRonize.hAMRonization_CONSOLIDATED_TERMS


        # AMR Term Consolidation
        File? amrtermconsolidation_isna = run_hAMRonize.amrtermconsolidation_isna
        File? amrtermconsolidation_all = run_hAMRonize.amrtermconsolidation_all
        File? amrtermconsolidation_over98 = run_hAMRonize.amrtermconsolidation_over98
        File? amrtermconsolidation_allidentity = run_hAMRonize.amrtermconsolidation_allidentity
    }

}

