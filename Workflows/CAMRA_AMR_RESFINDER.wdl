version 1.0

import "../Tasks/task_amrfinder.wdl" as amrfinder
import "../Tasks/task_abricate.wdl" as abricate
import "../Tasks/task_hamronization.wdl" as hamronize
import "../Tasks/task_resfinder.wdl" as resfinder
import "../Tasks/task_utilities.wdl" as utilities
import "../Tasks/task_rgi.wdl" as rgi
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
        String sample_name
        String organism

    }

    call resfinder.run_ResFinder{
        input:
            assembly = assembly,
            organism = organism

    }




    output {

        # ResFinder
        String resfinder_version = run_ResFinder.resfinder_version
        String resfinder_kma_version = run_ResFinder.resfinder_kma_version
        String resfinder_db_version = run_ResFinder.resfinder_db_version
        File resfider_asm_output = run_ResFinder.resfider_asm_output
        File resfinder_asm_hits = run_ResFinder.resfinder_asm_hits
        File resfinder_asm_argseq = run_ResFinder.resfinder_asm_argseq

    }

}

