version 1.0

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

     #BVBRC docker and build task to send assembly and conduct genome analysis & annotation
    call bvbrc.run_BVBRC_annotation_analysis {
        input:
            username = BVBRC_username,
            password = BVBRC_password,
            contigs_file_local = assembly,
            sample_name = sample_name,
            scientific_name = organism,  # "Genus species" from MASH, Optional
            taxonomy_id = 2
    }





    output {
       #BVBRC Annotation

        File bvbrc_annot_full_genome_report                 = run_BVBRC_annotation_analysis.bvbrc_annot_full_genome_report
        File bvbrc_annot_genome_annotation                  = run_BVBRC_annotation_analysis.bvbrc_annot_genome_annotation
        File bvbrc_annot_amr_annotation                     = run_BVBRC_annotation_analysis.bvbrc_annot_amr_annotation
        File bvbrc_annot_quality                            = run_BVBRC_annotation_analysis.bvbrc_annot_quality
        File bvbrc_annot_transformed_amrhits                = run_BVBRC_annotation_analysis.bvbrc_annot_transformed_amrhits
        File bvbrc_annot_transformed_predictedresistance    = run_BVBRC_annotation_analysis.bvbrc_annot_transformed_predictedresistance
        File bvbrc_annot_feature_protein                    = run_BVBRC_annotation_analysis.bvbrc_annot_feature_protein
    }

}
