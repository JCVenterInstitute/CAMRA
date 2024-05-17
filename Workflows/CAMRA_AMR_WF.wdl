version 1.0

import "../Tasks/AMR_finder.wdl" as amrfinder
import "../Tasks/Abricate.wdl" as abricate
import "../Tasks/hARMonization.wdl" as hamronize
import "../Tasks/ResFinder.wdl" as resfinder


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
        # File? pubmlst_DB
        File plasmidfinder_DB
        Array[File] AMR_files
        String plasmidfinder_DB
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

    call hamronize.run_Hamronize{
        input:
            sample_name = sample_name,
            #AMR_files = AMR_files
            AMR_files = [run_Abricate.abricate_ncbiDB_tsv_output, run_Abricate.abricate_cardDB_tsv_output,run_Abricate.abricate_resfinderDB_tsv_output, run_Abricate.abricate_vfdb_tsv_output, run_Abricate.abricate_argannotBD_tsv_output, run_AMRfinderPlus.AMRfinder_txt_output ]
    }
    
    call resfinder.run_ResFinder{
        input:
            assembly = assembly,
            read1 = read1,
            read2 = read2
        
    }
    
    
}
