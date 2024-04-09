version 1.0

import "../Tasks/Analysis/AMR_finder.wdl" as amrfinder
import "../Tasks/Analysis/MLST.wdl" as mlst


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
        File pubmlst_DB
    }

    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            sample_name = sample_name,
            organism = organism
    }
    call mlst.run_MLST{
        input:
            assembly = assembly,
            sample_name = sample_name,
            organism = organism,
            pubmlst_DB = pubmlst_DB
    }
    
}
