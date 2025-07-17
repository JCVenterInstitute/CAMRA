version 1.0

import "../Tasks/task_amrfinder.wdl" as amrfinder

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
         String sample_name
       }

    # Task to combine genus and species



    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            organism = organism
    }

   output {
         # AMR finder
        File amrfinder_all_output = run_AMRfinderPlus.amrfinder_all_output
        File amrfinder_stress_output = run_AMRfinderPlus.amrfinder_stress_output
        File amrfinder_virulence_output = run_AMRfinderPlus.amrfinder_virulence_output
        File amrfinder_amr_output = run_AMRfinderPlus.amrfinder_amr_output

        String amrinder_version = run_AMRfinderPlus.amrinder_version
        String amrfinder_db_version = run_AMRfinderPlus.amrfinder_db_version
        String amrfinder_date = run_AMRfinderPlus.amrfinder_date
    }

}
