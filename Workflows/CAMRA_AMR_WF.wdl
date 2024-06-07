version 1.0

import "../Tasks/AMR_finder.wdl" as amrfinder
import "../Tasks/Abricate.wdl" as abricate
import "../Tasks/hARMonization.wdl" as hamronize
import "../Tasks/ResFinder.wdl" as resfinder
import "../Tasks/AMR_harmonization.wdl" as amrharmonizer


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
        File read1
        File read2
        String sample_name
        String organism 

    }

    call amrfinder.run_AMRfinderPlus {
        input:
            assembly = assembly,
            sample_name = sample_name,
            organism = organism
    }
    
    call abricate.run_Abricate{
        input:
            assembly = assembly,
            sample_name = sample_name,
    }
    call resfinder.run_ResFinder{
        input:
            assembly = assembly,
            read1 = read1,
            read2 = read2,
            organism = organism
        
    }
    call hamronize.run_Hamronize{
        input:

            AMR_files = [run_Abricate.abricate_ncbiDB_tsv_output,
            run_Abricate.abricate_cardDB_tsv_output, 
            run_Abricate.abricate_resfinderDB_tsv_output, 
            run_Abricate.abricate_argannotBD_tsv_output, 
            run_AMRfinderPlus.amrfinder_amr_output,
            run_ResFinder.resfider_asm_output,
            run_ResFinder.resfinder_read_output],

            VIR_files = [run_Abricate.abricate_vfdb_tsv_output, 
            run_AMRfinderPlus.amrfinder_virulence_output]
    }
    
    call amrharmonizer.run_AMR_harmonizer {
        input:
        hamronize_amr_output = run_Hamronize.hAMRonization_amr_output
    }

    
    
}
