version 1.0 

task run_merqury {
    input {
        File assembly
        File read1
        File read2
        String sample_name
        String asm_size
        #TODO add option to change docker version
    }
    runtime{
        docker: 'miramastoras/merqury:latest'
        memory: "4G" #increasing memmory worked
    }
    command <<<
        #TODO add versioning
        # finding genome size and best kmer size

        total_length=~{asm_size} 
        best_k=$(best_k.sh $total_length) 
        best_k=$(echo "$best_k" | tail -n 1) 
        #Preparing meryl dbs
        meryl k=$best_k count output read1.meryl ~{read1}
        meryl k=$best_k count output read2.meryl ~{read2}

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
        Array[String] stdout_values = read_lines(stdout()) 
        String merqury_qv = stdout_values[11]
        String merqury_comp = stdout_values[12]
        File merqury_qv_file = "merqury_output/~{sample_name}.qv"
        File merqury_completeness_file = "merqury_output/~{sample_name}.completeness.stats"
    }
}