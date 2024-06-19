version 1.0

task run_fastQC {
    meta {
        # TODO
    }

    input {
        File read1
        File read2
    }

    runtime {
        docker: 'staphb/fastqc:0.12.1'
    }

    command <<<
        mkdir fastQC_output
        # echo $(fastqc --version 2>&1) | sed 's/abricate //' | tee VERSION

        # Run FastQC
        fastqc -o fastQC_output -j /usr/bin/java -f fastq ~{read1} ~{read2}

        # Define patterns
        pattern1="*R1_001_fastqc.zip"
        pattern2="*R2_001_fastqc.zip"

        # Find and unzip the files
        file1=$(find fastQC_output -name "$pattern1" -print -quit)
        file2=$(find fastQC_output -name "$pattern2" -print -quit)

        if [ -n "$file1" ]; then
            unzip "$file1" -d fastQC_output
        else
            echo "File matching pattern1 not found!"
        fi

        if [ -n "$file2" ]; then
            unzip "$file2" -d fastQC_output
        else
            echo "File matching pattern2 not found!"
        fi

        cp fastQC_output/*R1_001_fastqc/fastqc_data.txt fastQC_output/fastqc_R1_data.txt
        cp fastQC_output/*R1_001_fastqc/fastqc_report.html fastQC_output/fastqc_R1_report.html
        cp fastQC_output/*R1_001_fastqc/summary.txt fastQC_output/fastqc_R1_summary.txt

        cp fastQC_output/*R2_001_fastqc/fastqc_data.txt fastQC_output/fastqc_R2_data.txt
        cp fastQC_output/*R2_001_fastqc/fastqc_report.html fastQC_output/fastqc_R2_report.html
        cp fastQC_output/*R2_001_fastqc/summary.txt fastQC_output/fastqc_R2_summary.txt

        # Count the occurrences of 'PASS', 'WARN', and 'FAIL' in the R1 summary file
        pass_count_R1=$(grep -o 'PASS' fastQC_output/fastqc_R1_summary.txt | wc -l)
        warn_count_R1=$(grep -o 'WARN' fastQC_output/fastqc_R1_summary.txt | wc -l)
        fail_count_R1=$(grep -o 'FAIL' fastQC_output/fastqc_R1_summary.txt | wc -l)
        result_string_R1="P:${pass_count_R1} | W:${warn_count_R1} | F:${fail_count_R1}"
        echo $result_string_R1 > fastQC_output/result_string_R1.txt

        pass_count_R2=$(grep -o 'PASS' fastQC_output/fastqc_R2_summary.txt | wc -l)
        warn_count_R2=$(grep -o 'WARN' fastQC_output/fastqc_R2_summary.txt | wc -l)
        fail_count_R2=$(grep -o 'FAIL' fastQC_output/fastqc_R2_summary.txt | wc -l)
        result_string_R2="P:${pass_count_R2} | W:${warn_count_R2} | F:${fail_count_R2}"
        echo $result_string_R2 > fastQC_output/result_string_R2.txt

        # Extract values from *_data.txt
        total_sequences_R1=$(grep 'Total Sequences' fastQC_output/fastqc_R1_data.txt | cut -f2)
        total_bases_R1=$(grep 'Total Bases' fastQC_output/fastqc_R1_data.txt | cut -f2)
        poor_quality_R1=$(grep 'Sequences flagged as poor quality' fastQC_output/fastqc_R1_data.txt | cut -f2)
        sequence_length_R1=$(grep 'Sequence length' fastQC_output/fastqc_R1_data.txt | cut -f2)
        gc_content_R1=$(grep '%GC' fastQC_output/fastqc_R1_data.txt | cut -f2)

        total_sequences_R2=$(grep 'Total Sequences' fastQC_output/fastqc_R2_data.txt | cut -f2)
        total_bases_R2=$(grep 'Total Bases' fastQC_output/fastqc_R2_data.txt | cut -f2)
        poor_quality_R2=$(grep 'Sequences flagged as poor quality' fastQC_output/fastqc_R2_data.txt | cut -f2)
        sequence_length_R2=$(grep 'Sequence length' fastQC_output/fastqc_R2_data.txt | cut -f2)
        gc_content_R2=$(grep '%GC' fastQC_output/fastqc_R2_data.txt | cut -f2)

        # Save values to files
        echo $total_sequences_R1 > fastQC_output/total_sequences_R1.txt
        echo $total_bases_R1 > fastQC_output/total_bases_R1.txt
        echo $poor_quality_R1 > fastQC_output/poor_quality_R1.txt
        echo $sequence_length_R1 > fastQC_output/sequence_length_R1.txt
        echo $gc_content_R1 > fastQC_output/gc_content_R1.txt

        echo $total_sequences_R2 > fastQC_output/total_sequences_R2.txt
        echo $total_bases_R2 > fastQC_output/total_bases_R2.txt
        echo $poor_quality_R2 > fastQC_output/poor_quality_R2.txt
        echo $sequence_length_R2 > fastQC_output/sequence_length_R2.txt
        echo $gc_content_R2 > fastQC_output/gc_content_R2.txt



    >>>

    output {
        File fastQC_R1_data = "fastQC_output/fastqc_R1_data.txt"
        File fastQC_R2_data = "fastQC_output/fastqc_R2_data.txt"
        File fastQC_R1_html = "fastQC_output/fastqc_R1_report.html"
        File fastQC_R2_html = "fastQC_output/fastqc_R2_report.html"
        File fastQC_R1_summary = "fastQC_output/fastqc_R1_summary.txt"
        File fastQC_R2_summary = "fastQC_output/fastqc_R2_summary.txt"

        String fastQC_R1_PassWarnFail = read_string("fastQC_output/result_string_R1.txt")
        String fastQC_R2_PassWarnFail = read_string("fastQC_output/result_string_R2.txt")

        String fastQC_R1_total_sequences = read_string("fastQC_output/total_sequences_R1.txt")
        String fastQC_R1_total_bases = read_string("fastQC_output/total_bases_R1.txt")
        String fastQC_R1_poor_quality = read_string("fastQC_output/poor_quality_R1.txt")
        String fastQC_R1_sequence_length = read_string("fastQC_output/sequence_length_R1.txt")
        String fastQC_R1_gc_content = read_string("fastQC_output/gc_content_R1.txt")
        String fastQC_R2_total_sequences = read_string("fastQC_output/total_sequences_R2.txt")
        String fastQC_R2_total_bases = read_string("fastQC_output/total_bases_R2.txt")
        String fastQC_R2_poor_quality = read_string("fastQC_output/poor_quality_R2.txt")
        String fastQC_R2_sequence_length = read_string("fastQC_output/sequence_length_R2.txt")
        String fastQC_R2_gc_content = read_string("fastQC_output/gc_content_R2.txt")
        
        # String fastQC_version = read_string("VERSION")
    }
}

