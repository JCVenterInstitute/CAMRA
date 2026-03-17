version 1.0

import "../Tasks/task_bvbrc.wdl" as bvbrc

workflow amr_analysis_no_assembly   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Analysis on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
       String BVBRC_username
        String BVBRC_password
        String sample_name
        String date
    }
    call bvbrc.run_BVBRC_get_AMR_results {
        input: 
            username = BVBRC_username,
            password = BVBRC_password,
            date = date,
            sample_name = sample_name
    }
    output {
       #BVBRC Annotation

        File bvbrc_annot_genome_annotation                  = run_BVBRC_get_AMR_results.bvbrc_annot_genome_annotation
        File bvbrc_annot_amr_annotation                     = run_BVBRC_get_AMR_results.bvbrc_annot_amr_annotation
        File bvbrc_annot_quality                            = run_BVBRC_get_AMR_results.bvbrc_annot_quality
        File bvbrc_annot_transformed_amrhits                = run_BVBRC_get_AMR_results.bvbrc_annot_transformed_amrhits
        File bvbrc_annot_transformed_predictedresistance    = run_BVBRC_get_AMR_results.bvbrc_annot_transformed_predictedresistance
    }
}
