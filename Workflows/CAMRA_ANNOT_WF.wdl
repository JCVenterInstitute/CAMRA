version 1.0

import "../Tasks/plasmidfinder.wdl" as plasmidfinder

workflow assembly_analysis   {
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

}