version 1.0

import "../Tasks/BV-BRC_tasks.wdl" as bvbrc

workflow test_wdl  {
    input {
        File annotation_genome
        File genome_amr_json
        File quality_json

    }
    call bvbrc.run_hamronize_reformatting {
        input:
            annotation_genome = annotation_genome,
            genome_amr_json = genome_amr_json,
            quality_json = quality_json
    }

    output{
    # BVBRC output for harmonization 
        File bvbrc_amr_annotation = run_hamronize_reformatting.bvbrc_amr_annotation
        File bvbrc_predicted_resistance = bvbrc_predicted_resistance 
    }  
}
