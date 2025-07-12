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
        # File? read1
        # File? read2
        File? blast_query
        String sample_name
        String organism
        File? bvbrc_amr_file

        String? bvbrc_assembly_path
        String BVBRC_username
        String BVBRC_password
        String? bvbrc_timestamp
    }

    # Task to combine genus and species



    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            organism = organism
    }

    call rgi.run_RGI {
        input:
            assembly = assembly
    }

    if (defined(blast_query)){
        call amrfinder.run_Query_Blastn {
            input:
                assembly = assembly,
                query = blast_query
        }
    }

    call abricate.run_Abricate{
        input:
            assembly = assembly,
            sample_name = sample_name
    }

    call resfinder.run_ResFinder{
        input:
            assembly = assembly,
            organism = organism

    }

     #BVBRC docker and build task to send assembly and conduct genome analysis & annotation
    call bvbrc.run_BVBRC_annotation_analysis {
        input:
            username = BVBRC_username,
            password = BVBRC_password,
            bvbrc_assembly_path = bvbrc_assembly_path,
            contigs_file_local = assembly,
            sample_name = sample_name,
            timestamp = bvbrc_timestamp,
            scientific_name = organism,  # "Genus species" from MASH, Optional
            taxonomy_id = 2
    }


    call hamronize.run_hAMRonize {
        input:
            assembly = assembly,
            # AMR Output
            # abricate_ncbiDB_tsv_output = run_Abricate.abricate_ncbiDB_tsv_output,
            abricate_cardDB_tsv_output  = run_Abricate.abricate_cardDB_tsv_output,
            # abricate_resfinderDB_tsv_output  = run_Abricate.abricate_resfinderDB_tsv_output,
            abricate_argannotDB_tsv_output  = run_Abricate.abricate_argannotDB_tsv_output,

            amrfinder_amr_output = run_AMRfinderPlus.amrfinder_amr_output,

            resfider_asm_output = run_ResFinder.resfider_asm_output,

            rgi_CARD_diamond_tsv_output = run_RGI.rgi_CARD_diamond_tsv_output,
            rgi_CARD_blast_tsv_output = run_RGI.rgi_CARD_blast_tsv_output,

            #bvbrc_amr_file = run_BVBRC_annotation_analysis.bvbrc_annot_transformed_amrhits,

            # Virulence Output

            VIR_files = [
                run_Abricate.abricate_vfdb_tsv_output,
                run_AMRfinderPlus.amrfinder_virulence_output
            ]
        }



    output {
        # BVBRC Annotation

       # File bvbrc_annot_full_genome_report                 = run_BVBRC_annotation_analysis.bvbrc_annot_full_genome_report
       # File bvbrc_annot_genome_annotation                  = run_BVBRC_annotation_analysis.bvbrc_annot_genome_annotation
       # File bvbrc_annot_amr_annotation                     = run_BVBRC_annotation_analysis.bvbrc_annot_amr_annotation
       # File bvbrc_annot_quality                            = run_BVBRC_annotation_analysis.bvbrc_annot_quality
       # File bvbrc_annot_transformed_amrhits                = run_BVBRC_annotation_analysis.bvbrc_annot_transformed_amrhits
       # File bvbrc_annot_transformed_predictedresistance    = run_BVBRC_annotation_analysis.bvbrc_annot_transformed_predictedresistance
       # File bvbrc_annot_feature_protein                    = run_BVBRC_annotation_analysis.bvbrc_annot_feature_protein

        # Optional Output - blast against userinput query
        File? blastn_output = run_Query_Blastn.blastn_output

        # AMR finder
        File amrfinder_all_output = run_AMRfinderPlus.amrfinder_all_output
        File amrfinder_stress_output = run_AMRfinderPlus.amrfinder_stress_output
        File amrfinder_virulence_output = run_AMRfinderPlus.amrfinder_virulence_output
        File amrfinder_amr_output = run_AMRfinderPlus.amrfinder_amr_output

        String amrinder_version = run_AMRfinderPlus.amrinder_version
        String amrfinder_db_version = run_AMRfinderPlus.amrfinder_db_version
        String amrfinder_date = run_AMRfinderPlus.amrfinder_date

        # Abricate
        # File abricate_ncbiDB_tsv_output = run_Abricate.abricate_ncbiDB_tsv_output
        File abricate_cardDB_tsv_output = run_Abricate.abricate_cardDB_tsv_output
        # File abricate_resfinderDB_tsv_output = run_Abricate.abricate_resfinderDB_tsv_output
        File abricate_vfdb_tsv_output = run_Abricate.abricate_vfdb_tsv_output
        File abricate_argannotDB_tsv_output = run_Abricate.abricate_argannotDB_tsv_output

        File abricate_DB_version = run_Abricate.abricate_DB_version
        String abricate_version = run_Abricate.abricate_version
        String abricate_date = run_Abricate.abricate_date

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

        # ResFinder
        String resfinder_version = run_ResFinder.resfinder_version
        String resfinder_kma_version = run_ResFinder.resfinder_kma_version
        String resfinder_db_version = run_ResFinder.resfinder_db_version
        File resfider_asm_output = run_ResFinder.resfider_asm_output
        File resfinder_asm_hits = run_ResFinder.resfinder_asm_hits
        File resfinder_asm_argseq = run_ResFinder.resfinder_asm_argseq


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

