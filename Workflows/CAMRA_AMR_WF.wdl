version 1.0

import "../Tasks/task_amrfinder.wdl" as amrfinder
import "../Tasks/task_abricate.wdl" as abricate
import "../Tasks/task_hamronization.wdl" as hamronize
import "../Tasks/task_resfinder.wdl" as resfinder
import "../Tasks/task_utilities.wdl" as utilities
import "../Tasks/task_rgi.wdl" as rgi


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
        File read1
        File read2
        File blast_carbapen_query
        String sample_name
        String genus
        String species
    }

    # Task to combine genus and species
    call utilities.run_taxajoin {
        input:
            genus = genus,
            species = species
    }


    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            organism = run_taxajoin.organism 
    }

    call rgi.run_RGI {
        input:
            assembly = assembly
    }

    call amrfinder.run_Query_Blastn {
        input:
            assembly = assembly,
            query = blast_carbapen_query
    }
    
    call abricate.run_Abricate{
        input:
            assembly = assembly,
            sample_name = sample_name,
    }
    call resfinder.run_ResFinder{
        input:
            assembly = assembly,
            read1 = read1,
            read2 = read2,
            organism = run_taxajoin.organism
        
    }
    call hamronize.run_hAMRonize{
        input:

            AMR_files = [
                run_Abricate.abricate_ncbiDB_tsv_output,
                run_Abricate.abricate_cardDB_tsv_output, 
                run_Abricate.abricate_resfinderDB_tsv_output, 
                run_Abricate.abricate_argannotBD_tsv_output, 

                run_AMRfinderPlus.amrfinder_amr_output,

                run_ResFinder.resfider_asm_output,
                run_ResFinder.resfinder_read_output,

                run_RGI.rgi_CARD_diamond_tsv_output,
                run_RGI.rgi_CARD_blast_tsv_output
            ],

            VIR_files = [
                run_Abricate.abricate_vfdb_tsv_output, 
                run_AMRfinderPlus.amrfinder_virulence_output
            ]
    }
    


    output{
        # Optional Output - blast against userinput query
        File blastn_output = run_Query_Blastn.blastn_output

        # AMR finder
        File amrfinder_all_output = run_AMRfinderPlus.amrfinder_all_output
        File amrfinder_stress_output = run_AMRfinderPlus.amrfinder_stress_output
        File amrfinder_virulence_output = run_AMRfinderPlus.amrfinder_virulence_output
        File amrfinder_amr_output = run_AMRfinderPlus.amrfinder_amr_output

        String amrinder_version = run_AMRfinderPlus.amrinder_version
        String amrfinder_db_version = run_AMRfinderPlus.amrfinder_db_version
        String amrfinder_date = run_AMRfinderPlus.amrfinder_date

        # Abricate 
        File abricate_ncbiDB_tsv_output = run_Abricate.abricate_ncbiDB_tsv_output
        File abricate_cardDB_tsv_output = run_Abricate.abricate_cardDB_tsv_output
        File abricate_resfinderDB_tsv_output = run_Abricate.abricate_resfinderDB_tsv_output
        File abricate_vfdb_tsv_output = run_Abricate.abricate_vfdb_tsv_output 
        File abricate_argannotBD_tsv_output = run_Abricate.abricate_argannotBD_tsv_output

        File abricate_DB_version = run_Abricate.abricate_DB_version
        String abricate_version = run_Abricate.abricate_version 
        String abricate_date = run_Abricate.abricate_date 

        # hAMRonization
        String hAMRonization_version = run_hAMRonize.hAMRonization_version
        String hAMRonization_date = run_hAMRonize.hAMRonization_date
        File hAMRonization_amr_output = run_hAMRonize.hAMRonization_amr_output
        File hAMRonization_vir_output = run_hAMRonize.hAMRonization_vir_output

        # AMR Term Consolidation
        File amrtermconsolidation_isna = run_hAMRonize.amrtermconsolidation_isna
        File amrtermconsolidation_all = run_hAMRonize.amrtermconsolidation_all
        File amrtermconsolidation_over98 = run_hAMRonize.amrtermconsolidation_over98 
        File amrtermconsolidation_allidentity = run_hAMRonize.amrtermconsolidation_allidentity

        # ResFinder
        String resfinder_version = run_ResFinder.resfinder_version
        String resfinder_kma_version = run_ResFinder.resfinder_kma_version
        String resfinder_db_version = run_ResFinder.resfinder_db_version
        File resfider_asm_output = run_ResFinder.resfider_asm_output
        File resfinder_read_output = run_ResFinder.resfinder_read_output
        File resfinder_asm_hits = run_ResFinder.resfinder_asm_hits
        File resfinder_read_hits = run_ResFinder.resfinder_read_hits
        File resfinder_asm_argseq = run_ResFinder.resfinder_asm_argseq
        File resfinder_read_argseq = run_ResFinder.resfinder_read_argseq

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
