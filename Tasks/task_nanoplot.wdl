version 1.0

task run_nanoplot {
  input {
    File ont_reads # intended for ONT data only
    String sample_name
  }
  command <<<
    # get date
    date | tee "DATE"
    # get version
    NanoPlot --version | tee "VERSION"
    
    mkdir ~{sample_name}

    NanoPlot \
        --fastq ~{ont_reads} \
        -o ~{sample_name}\
        -t 12 \
        --N50 \

    # grep read statistics from tsv stats file
    grep "Number of reads:" ~{sample_name}/NanoStats.txt | awk '{print $NF}' | tee NUMBER_OF_READS
    NUM_BASES=$(grep "Total bases:" ~{sample_name}/NanoStats.txt | awk '{print $NF}' | tee NUMBER_OF_BASES)
    grep "Median read length:" ~{sample_name}/NanoStats.txt | awk '{print $NF}'| tee MEDIAN_READ_LENGTH
    grep "Mean read length:" ~{sample_name}/NanoStats.txt | awk '{print $NF}' | tee MEAN_READ_LENGTH
    # grep "read_length_stdev" ~{sample_name}/NanoStats.txt | awk '{print $NF}' | tee READ_LENGTH_STDEV
    grep "Read length N50:" ~{sample_name}/NanoStats.txt | awk '{print $NF}' | tee N50
    grep "Mean read quality:" ~{sample_name}/NanoStats.txt | awk '{print $NF}' | tee MEAN_QUAL
    grep "Median read quality:" ~{sample_name}/NanoStats.txt | awk '{print $NF}'| tee MEDIAN_QUAL


  >>>
  output {
    File nanoplot_html = "~{sample_name}/NanoPlot-report.html"
    File nanoplot_tsv = "~{sample_name}/NanoStats.txt"
    String nanoplot_num_reads = read_string("NUMBER_OF_READS")
    String nanoplot_median_readlength = read_string("MEDIAN_READ_LENGTH")
    String nanoplot_mean_readlength = read_string("MEAN_READ_LENGTH")
    # Float nanoplot_stdev_readlength = read_float("READ_LENGTH_STDEV")
    String nanoplot_n50 = read_string("N50")
    String nanoplot_mean_quality = read_string("MEAN_QUAL")
    String nanoplot_median_quality = read_string("MEDIAN_QUAL")
    String nanoplot_version = read_string("VERSION")
    String nanoplot_date = read_string("DATE")
  }
  runtime {
    docker: "nanozoo/nanoplot:latest"
  }
}