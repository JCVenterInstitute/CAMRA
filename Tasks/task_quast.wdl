version 1.0

task run_Quast {
    meta {
        description: "Evaluates genome/metagenome assemblies by computing various metrics."
        gitrepository: "https://github.com/ablab/quast"
    }

    input {
        File assembly
        Int min_contigs = 500
    }

    runtime {
        docker: 'staphb/quast:5.2.0-slim'
    }

    command <<<
        date | tee DATE
        quast.py --version | tee VERSION
        quast.py -o quast_output ~{if defined(min_contigs) then "--min-contig " + min_contigs else ""} ~{assembly}
        # Extract specific metrics from the TSV report
        awk -F"\t" '/Largest contig/ {print $2}' quast_output/report.tsv > quast_output/largest_contig.txt
        awk -F'\t' '/^Total length\t/ {print $2}' quast_output/report.tsv > quast_output/total_length.txt
        awk -F"\t" '/N50/ {print $2}' quast_output/report.tsv  > quast_output/N50.txt
        awk -F"\t" '/N90/ {print $2}' quast_output/report.tsv > quast_output/N90.txt
        awk -F"\t" '/L50/ {print $2}' quast_output/report.tsv > quast_output/L50.txt
        awk -F"\t" '/L90/ {print $2}' quast_output/report.tsv > quast_output/L90.txt
    >>>

    output {
        String quast_version = read_string("VERSION")
        String quast_date = read_string("DATE")
        
        File quast_report = "quast_output/report.tsv"
        Int quast_contig_largest = read_int("quast_output/largest_contig.txt")
        Int quast_total_length = read_int("quast_output/total_length.txt")
        Int quast_N50 = read_int("quast_output/N50.txt")
        Int quast_N90 = read_int("quast_output/N90.txt")
        Int quast_L50 = read_int("quast_output/L50.txt")
        Int quast_L90 = read_int("quast_output/L90.txt")
    }
}