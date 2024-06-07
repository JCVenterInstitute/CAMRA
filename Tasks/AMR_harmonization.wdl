version 1.0


task run_AMR_harmonizer {
    input {
            File hamronize_amr_output
        }
    runtime {
        docker: 'danylmb/amrharmonize:1.0'
    }
    command <<<

        python /usr/src/app/amrharmonization.py ~{hamronize_amr_output}

        
    >>>

    output {
        File amrhamronization_isna = "hamronization_isna.tsv"
        File amrhamronization_all = "hamronization_all.tsv"
        File amrharmonization_over98 = "harmonized_amr_over98identity.tsv"
        File amrharmonization_all = "harmonized_amr_allidentity.tsv"
    }




}