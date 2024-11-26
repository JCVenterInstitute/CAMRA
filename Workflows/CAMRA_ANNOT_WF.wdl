version 1.0

import "../Tasks/plasmidfinder.wdl" as plasmidfinder
import "../Tasks/BV-BRC_tasks.wdl" as bvbrc

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
        File local_assembly_path
        String? bvbrc_assembly_path
        String sample_name
        String BVBRC_username
        String BVBRC_password
        String? timestamp
        String scientific_name
        File plasmidfinder_DB
    }

    call plasmidfinder.run_PlasmidFinder {
        #there is also plasflow and plasmidspades
        input:
            assembly = local_assembly_path,
            sample_name = sample_name,
            database = plasmidfinder_DB 
    }
    # BVBRC docker and build task to send assembly and conduct genome analysis & annotation
    call bvbrc.run_annotation_analysis {
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

    #TODO Add phagefinder, I am waiting to have the list of important output files. 


    # call pgap, prokka, bakta

    # FCS
    # IS elements
    # Mobile element Finder


    output {
        File full_genome_report = run_annotation_analysis.full_genome_report
        File annotation_genome_report = run_annotation_analysis.annotation_genome_report
        File annotation_xls = run_annotation_analysis.annotation_xls
        String plasmidfinder_plasmids_list = run_PlasmidFinder.plasmidfinder_plasmids_list
        String plasmidfinder_qty_hits =run_PlasmidFinder.plasmidfinder_qty_hits 
        File plasmidfinder_tsv_output = run_PlasmidFinder.plasmidfinder_tsv_output
        File plasmidfinder_seq_output = run_PlasmidFinder.plasmidfinder_seq_output
    }

}