version 1.0 

task run_merqury {
    input {
        File assembly
        File read1
        File read2
        Int asm_size
        String memory = "4G" 
    }

    runtime{
        docker : 'danylmb/merqury:1.4.1-build2'
        cpu: 4
        memory: "32 GB"
        disk: "/tmp 1000 SSD"
    }
    command <<<
        date | tee DATE
        echo "V1.3" | tee MERQURY_VERSION

        total_length=~{asm_size} 
        best_k=$(best_k.sh $total_length) 
        best_k=$(echo "$best_k" | tail -n 1) 
        #Preparing meryl dbs
        meryl k=$best_k count ~{'memory=' + memory} threads=1 output read1.meryl ~{read1}
        meryl k=$best_k count ~{'memory=' + memory} threads=1 output read2.meryl ~{read2}

        meryl union-sum output sample.meryl read1.meryl read2.meryl 

        #Using Merqury
        mkdir merqury_output && cd merqury_output
        merqury.sh ../sample.meryl ~{assembly} sample && echo "MERQURY DONE"
        #Outputs of Intrest
        qv_line=$(tail -n 3 sample.qv | head -n 1)
        qv_number=$(echo "$qv_line" | awk '{print $4}')
        echo $qv_number
        comp_line=$(tail -n 6 sample.completeness.stats | head -n 1)
        comp_number=$(echo "$comp_line" | awk '{print $5}')
        echo $comp_number

        >>>
    output {
        String merqury_version = read_string("MERQURY_VERSION")
        String merqury_date = read_string("DATE")

        Array[String] stdout_values = read_lines(stdout()) 

        String merqury_qv = stdout_values[13]
        String merqury_comp = stdout_values[14]

        File merqury_qv_file = "merqury_output/sample.qv"
        File merqury_completeness_file = "merqury_output/sample.completeness.stats"
    }
}