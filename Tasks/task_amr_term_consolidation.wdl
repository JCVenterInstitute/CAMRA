version 1.0


task run_AMR_Term_Consolidation {
    input {
            File hamronize_amr_output
        }
    runtime {
        docker: 'danylmb/amrharmonize:1.0-build1'
    }
    command <<<

        python /usr/src/app/amrharmonization.py ~{hamronize_amr_output}

        
    >>>

    output {
        File amrtermconsolidation_isna = "hamronization_isna.tsv"
        File amrtermconsolidation_all = "hamronization_all.tsv"
        File amrtermconsolidation_over98 = "harmonized_amr_over98identity.tsv"
        File amrtermconsolidation_allidentity = "harmonized_amr_allidentity.tsv"
    }




}