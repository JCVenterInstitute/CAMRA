version 1.0

import "../Tasks/task_mob_suite.wdl" as mob_suite

workflow run_mob_suite   {
    meta {
        author: "THomas Clarke"
        email: "tclarke@jcvi.org"
        description: "Run Mob Suite on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        File assembly
    }

    # Task to combine genus and species



    call mob_suite.run_mob_suite {
        input:
            assembly = assembly
    }
    
    output {
        File mob_contig_out = run_mob_suite.mob_contig_out
        File mob_typer_out = run_mob_suite.mob_typer_out
        File mob_typing_out = run_mob_suite.mob_typing_out
        }
}
