version 1.0

import "../Tasks/task_plasmidfinder.wdl" as plasmidfinder
import "../Tasks/task_bvbrc.wdl" as bvbrc
import "../Tasks/task_pgap.wdl" as pgap

workflow annotation_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Annotation on an Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        String sample_name
        String BVBRC_username
        String BVBRC_password

        String scientific_name # Example "Escherichia coli"

        # Databases
        File plasmidfinder_DB
        File pgap_supplemental_data

        # Use if the assembly is not stored in BV-BRC (i.e. locally or cloud)
        File local_assembly_path

        # Use if the assembly is stored in BV-BRC
        String? timestamp
        String? bvbrc_assembly_path
    }

    call plasmidfinder.run_PlasmidFinder {
        #there is also plasflow and plasmidspades
        input:
            assembly = local_assembly_path,
            sample_name = sample_name,
            database = plasmidfinder_DB 
    }

    # BVBRC docker and build task to send assembly and conduct genome analysis & annotation
    call bvbrc.run_BVBRC_annotation_analysis {
        input:
            username = BVBRC_username,
            password = BVBRC_password,
            bvbrc_assembly_path = bvbrc_assembly_path,
            contigs_file_local = local_assembly_path,
            sample_name = sample_name,
            timestamp = timestamp,
            scientific_name = scientific_name,  # "Genus species" from MASH, Optional
            taxonomy_id = 2
    }    

    call pgap.pgap_annotate {
    input:
        sample_id = sample_name,
        assembly_fasta = local_assembly_path,
        supplemental_data = pgap_supplemental_data, # Refer to https://github.com/ncbi/pgap/wiki/Installation
        taxon = scientific_name, # Example "Escherichia coli"
        pgap_version = "2024-07-18.build7555"
    }

    #TODO Add phagefinder, I am waiting to have the list of important output files. 



    # TODO - i might just have do install a MGE insetad of IS as MGE is the overcumposing umbrella
    # IS elements - Small, autonomous transposable elements (~700–2500 bp) with terminal inverted repeats and a transposase gene.

    # Mobile Genetic Elements (MGEs) - Broad category including IS, transposons, plasmids, phages, and integrative elements.
        # Find the diffrence between :
        #1. MGEFinder: Detects putative mobile genetic elements from sequencing reads.
            # 115 github start 
            # pubmed: 80 citated
        #2. MobileElementFinder: Identifies MGEs such as transposons and insertion sequences in bacterial genomes.
            # pubmed: 224 cited
        # Run both tools and compare their outputs. Pick the one that seems the best/ most informative. Check their use in papers. 


    output {

        File bvbrc_annot_full_genome_report = run_BVBRC_annotation_analysis.bvbrc_annot_full_genome_report 
        File bvbrc_annot_genome_annotation = run_BVBRC_annotation_analysis.bvbrc_annot_genome_annotation
        File bvbrc_annot_amr_annotation  = run_BVBRC_annotation_analysis.bvbrc_annot_amr_annotation
        File bvbrc_annot_quality = run_BVBRC_annotation_analysis.bvbrc_annot_quality
        File bvbrc_annot_transformed_amrhits =  run_BVBRC_annotation_analysis.bvbrc_annot_transformed_amrhits
        File bvbrc_annot_transformed_predictedresistance = run_BVBRC_annotation_analysis.bvbrc_annot_transformed_predictedresistance
        File bvbrc_annot_feature_protein = run_BVBRC_annotation_analysis.bvbrc_annot_feature_protein

        File? pgap_ani_tax_report = pgap_annotate.pgap_ani_tax_report
        File? pgap_ani_tax_report_xml = pgap_annotate.pgap_ani_tax_report_xml
        File? pgap_annotated_proteins_fasta = pgap_annotate.pgap_annotated_proteins_fasta
        File? pgap_genomic_fasta = pgap_annotate.pgap_genomic_fasta
        File? pgap_gbk_annotation = pgap_annotate.pgap_gbk_annotation
        File? pgap_gff3_annotation = pgap_annotate.pgap_gff3_annotation
        File? pgap_annotated_sequence = pgap_annotate.pgap_annotated_sequence
        File? pgap_annotated_cds = pgap_annotate.pgap_annotated_cds
        File? pgap_translated_cds = pgap_annotate.pgap_translated_cds
        File? pgap_gff_with_sequence = pgap_annotate.pgap_gff_with_sequence
        File? pgap_foreign_sequence = pgap_annotate.pgap_foreign_sequence
        File? pgap_completeness_contamination = pgap_annotate.pgap_completeness_contamination
        File? pgap_log = pgap_annotate.pgap_log

        String plasmidfinder_plasmids_list = run_PlasmidFinder.plasmidfinder_plasmids_list
        String plasmidfinder_qty_hits =run_PlasmidFinder.plasmidfinder_qty_hits 
        File plasmidfinder_tsv_output = run_PlasmidFinder.plasmidfinder_tsv_output
        File plasmidfinder_seq_output = run_PlasmidFinder.plasmidfinder_seq_output
    }

}