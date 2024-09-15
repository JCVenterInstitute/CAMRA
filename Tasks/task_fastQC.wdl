version 1.0

task run_fastQC {
    meta {
    description: "FastQC aims to provide a simple way to do some quality control checks on raw read data coming from high throughput sequencing pipelines."
    gitrepository: "https://github.com/s-andrews/FastQC"
    docker:"https://hub.docker.com/r/staphb/fastqc"
    cite:"Andrews, S. (2010), 'FASTQC. A quality control tool for high throughput read data' ."
    }

    input {
        File read1
        File read2
    }

    runtime {
        docker: 'staphb/fastqc:0.12.1'
        memory: "4 GB"
        cpu: 8
    }

    command <<<
        mkdir read1 read2
        echo $(fastqc --version) | sed 's/FastQC //' | tee VERSION
        date | tee DATE

        # Run FastQC
        fastqc -t 8 -o read1 -j /usr/bin/java -f fastq ~{read1} 
        fastqc -t 8 -o read2 -j /usr/bin/java -f fastq ~{read2}

        mv read1/*.html read1/read1.html
        mv read1/*.zip read1/read1.zip
        unzip -d read1 read1/read1.zip

        directories=$(find read1 -maxdepth 1 -type d | sed 1d)
        if [ $(echo "$directories" | wc -l) -eq 1 ]; then
            mv $directories read1_output
        else
            echo "More than one directory found"
        fi

        mv read2/*.html read2/read2.html
        mv read2/*.zip read2/read2.zip
        unzip -d read2 read2/read2.zip

        directories=$(find read2 -maxdepth 1 -type d | sed 1d)
        if [ $(echo "$directories" | wc -l) -eq 1 ]; then
            mv $directories read2_output
        else
            echo "More than one directory found"
        fi

        # Count the occurrences of 'PASS', 'WARN', and 'FAIL' in the R1 summary file
        pass_count_R1=$(grep -o 'PASS' read1_output/summary.txt | wc -l)
        warn_count_R1=$(grep -o 'WARN' read1_output/summary.txt | wc -l)
        fail_count_R1=$(grep -o 'FAIL' read1_output/summary.txt | wc -l)
        result_string_R1="P:${pass_count_R1} | W:${warn_count_R1} | F:${fail_count_R1}"
        echo $result_string_R1 > read1_output/result_string_R1.txt

        pass_count_R2=$(grep -o 'PASS' read2_output/summary.txt | wc -l)
        warn_count_R2=$(grep -o 'WARN' read2_output/summary.txt | wc -l)
        fail_count_R2=$(grep -o 'FAIL' read2_output/summary.txt | wc -l)
        result_string_R2="P:${pass_count_R2} | W:${warn_count_R2} | F:${fail_count_R2}"
        echo $result_string_R2 > read2_output/result_string_R2.txt

        # Extract values from *_data.txt
        total_reads_R1=$(grep 'Total Sequences' read1_output/fastqc_data.txt | cut -f2)
        total_bases_R1=$(grep 'Total Bases' read1_output/fastqc_data.txt | cut -f2)
        poor_quality_R1=$(grep 'Sequences flagged as poor quality' read1_output/fastqc_data.txt | cut -f2)
        read_length_R1=$(grep 'Sequence length' read1_output/fastqc_data.txt | cut -f2)
        gc_content_R1=$(grep '%GC' read1_output/fastqc_data.txt | cut -f2)

        total_reads_R2=$(grep 'Total Sequences' read2_output/fastqc_data.txt | cut -f2)
        total_bases_R2=$(grep 'Total Bases' read2_output/fastqc_data.txt | cut -f2)
        poor_quality_R2=$(grep 'Sequences flagged as poor quality' read2_output/fastqc_data.txt | cut -f2)
        read_length_R2=$(grep 'Sequence length' read2_output/fastqc_data.txt | cut -f2)
        gc_content_R2=$(grep '%GC' read2_output/fastqc_data.txt | cut -f2)

        # Save values to files
        echo $total_reads_R1 > read1_output/total_reads_R1.txt
        echo $total_bases_R1 > read1_output/total_bases_R1.txt
        echo $poor_quality_R1 > read1_output/poor_quality_R1.txt
        echo $read_length_R1 > read1_output/read_length_R1.txt
        echo $gc_content_R1 > read1_output/gc_content_R1.txt

        echo $total_reads_R2 > read2_output/total_reads_R2.txt
        echo $total_bases_R2 > read2_output/total_bases_R2.txt
        echo $poor_quality_R2 > read2_output/poor_quality_R2.txt
        echo $read_length_R2 > read2_output/read_length_R2.txt
        echo $gc_content_R2 > read2_output/gc_content_R2.txt



    >>>

    output {
        File fastQC_R1_html = "read1_output/fastqc_report.html"
        File fastQC_R2_html = "read2_output/fastqc_report.html"

        String fastQC_R1_PassWarnFail = read_string("read1_output/result_string_R1.txt")
        String fastQC_R2_PassWarnFail = read_string("read2_output/result_string_R2.txt")

        String fastQC_R1_total_reads = read_string("read1_output/total_reads_R1.txt")
        String fastQC_R1_total_bases = read_string("read1_output/total_bases_R1.txt")
        String fastQC_R1_poor_quality = read_string("read1_output/poor_quality_R1.txt")
        String fastQC_R1_read_length = read_string("read1_output/read_length_R1.txt")
        String fastQC_R1_gc_content = read_string("read1_output/gc_content_R1.txt")
        String fastQC_R2_total_reads = read_string("read2_output/total_reads_R2.txt")
        String fastQC_R2_total_bases = read_string("read2_output/total_bases_R2.txt")
        String fastQC_R2_poor_quality = read_string("read2_output/poor_quality_R2.txt")
        String fastQC_R2_read_length = read_string("read2_output/read_length_R2.txt")
        String fastQC_R2_gc_content = read_string("read2_output/gc_content_R2.txt")
        
        String fastQC_version = read_string("VERSION")
        String fastQC_date = read_string("DATE")
    }
}

