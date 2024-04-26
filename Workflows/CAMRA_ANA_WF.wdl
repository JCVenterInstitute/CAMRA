version 1.0

import "../Tasks/Analysis/AMR_finder.wdl" as amrfinder
import "../Tasks/Analysis/MLST.wdl" as mlst
import "../Tasks/Analysis/plasmidfinder.wdl" as plasmidfinder
import "../Tasks/Analysis/Abricate.wdl" as abricate
import "../Tasks/Analysis/hARMonization.wdl" as hamronize


workflow assembly_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Analysis on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        String sample_name
        File assembly
        String organism 
        # File? pubmlst_DB
        File plasmidfinder_DB
        #String plasmidfinder_DB
    }

    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            sample_name = sample_name,
            organism = organism
    }
    

    # call mlst.run_MLST {
    #     input:
    #         assembly = assembly,
    #         sample_name = sample_name
    #         #pubmlst_DB = pubmlst_DB
    # }

    # call plasmidfinder.run_PlasmidFinder {
    #     input:
    #         assembly = assembly,
    #         sample_name = sample_name,
    #         database = plasmidfinder_DB 
    # }

    call abricate.run_Abricate{
        input:
            assembly = assembly,
            sample_name = sample_name,
    }

    call hamronize.run_Hamronize{
        input:
            sample_name = sample_name,
            AMR_files = [run_Abricate.abricate_ncbiDB_tsv_output, run_Abricate.abricate_cardDB_tsv_output, run_AMRfinderPlus.AMRfinder_tsv_output ]
    }

    
    
}
