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




}