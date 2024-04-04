version 1.0

import "../Tasks/Analysis/AMR_finder.wdl" as amrfinder


workflow assembly_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Analysis on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genomes"
    }
    input {
        String sample_name
        # File read1
        # File read2
        File assembly
        String organism 
        #Q? where should I put in the docker image i want to default used if an option is not provided?
        String amrfinder_docker
    }

    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            sample_name = sample_name,
            organism = organism,
            docker = amrfinder_docker 
    }
}
