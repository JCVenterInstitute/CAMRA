version 1.0

import "../Tasks/task_utilities.wdl" as utilities



workflow amr_analysis   {
    meta {
        author: "Thomas Clarke"
        email: "tclarke@jcvi.org"
        description: "Run Analysis on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        String genus
        String species
    }

    # Task to combine genus and species
    call utilities.run_taxajoin {
        input:
            genus = genus,
            species = species
    }
    
    output {
    	String Organism = run_taxajoin.organism
    }
}
