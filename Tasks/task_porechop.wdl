version 1.0

task run_porechop {
  input {
    File ont_reads
    String sample_name
  }

  runtime {
      docker: "staphb/porechop:latest"
  }

  command <<<
    # get date
    date | tee DATE
    # get version   
    porechop --version | tee VERSION

    #run porechop
    porechop -i ~{ont_reads}   -o ~{sample_name}.fastq.gz      -t 12
  >>>

  output {
    File porechop_reads = "~{sample_name}.fastq.gz"
    String porechop_version = read_string("VERSION")
    String porechop_date = read_string("DATE")
  }
}