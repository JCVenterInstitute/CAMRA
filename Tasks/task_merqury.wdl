version 1.0 

task run_merqury {
    input {
        File assembly
        File read1
        File read2
        String sample_name
        String asm_size
        String? memory = "4G" 
    }

    runtime{
        docker : 'danylmb/merqury:1.4.1'
        memory: "4G" 
    }
    command <<<
        echo "V1.3" | tee MERQURY_VERSION

        total_length=~{asm_size} 
        best_k=$(best_k.sh $total_length) 
        best_k=$(echo "$best_k" | tail -n 1) 
        #Preparing meryl dbs
        meryl k=$best_k count ~{'memory=' + memory} threads=1 output read1.meryl ~{read1}
        meryl k=$best_k count ~{'memory=' + memory} threads=1 output read2.meryl ~{read2}

        meryl union-sum output ~{sample_name}.meryl read1.meryl read2.meryl 

        #Using Merqury
        mkdir merqury_output && cd merqury_output
        merqury.sh ../~{sample_name}.meryl ~{assembly} ~{sample_name} && echo "MERQURY DONE"
        #Outputs of Intrest
        qv_line=$(tail -n 3 ~{sample_name}.qv | head -n 1)
        qv_number=$(echo "$qv_line" | awk '{print $4}')
        echo $qv_number
        comp_line=$(tail -n 6 ~{sample_name}.completeness.stats | head -n 1)
        comp_number=$(echo "$comp_line" | awk '{print $5}')
        echo $comp_number

        >>>
    output {
        String merqury_version = read_string("MERQURY_VERSION")

        Array[String] stdout_values = read_lines(stdout()) 

        String merqury_qv = stdout_values[12]
        String merqury_comp = stdout_values[13]

        File merqury_qv_file = "merqury_output/~{sample_name}.qv"
        File merqury_completeness_file = "merqury_output/~{sample_name}.completeness.stats"
    }
}