version 1.0

import "../Tasks/plasmidfinder.wdl" as plasmidfinder

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
        File assembly
        File plasmidfinder_DB
    }

    call plasmidfinder.run_PlasmidFinder {
        input:
            assembly = assembly,
            sample_name = sample_name,
            database = plasmidfinder_DB 
    }

    output {

        String plasmidfinder_plasmids_list = run_PlasmidFinder.plasmidfinder_plasmids_list
        String plasmidfinder_qty_hits =run_PlasmidFinder.plasmidfinder_qty_hits 
        File plasmidfinder_tsv_output = run_PlasmidFinder.plasmidfinder_tsv_output
        File plasmidfinder_seq_output = run_PlasmidFinder.plasmidfinder_seq_output
    }

}