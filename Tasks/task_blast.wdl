
version 1.0

task blastn {
input {
    String samplename
    File? query
    File subject
    Float evalue = 0.001
    Int max_hsps = 1
    Int percent_identity = 90
    String docker = "staphb/blast:2.15.0"
}

command <<<
    # version
    tblastn -version | grep tblastn | cut -d " " -f2 > VERSION

    # tblastn execution
    hits=""
    if [[ -n "~{query}" && -n "~{subject}" ]]; then
        gunzip -c ~{query} > assembly.fasta

        blastn -query assembly.fasta -subject ~{subject} -evalue ~{evalue} -max_hsps ~{max_hsps} -out blast.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

        awk -v id="~{percent_identity}" 'BEGIN{print "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen"}; $3>=id {print}' blast.txt > ~{samplename}.tblastn.tsv

        hits=$(awk -v id="~{percent_identity}" '$3>=id {print $1}' blast.txt | tr '\n' ', ' | sed 's/,$//')

        rm assembly.fasta
    fi
    
    if [ -z "$hits" ]; then
        echo "No blast hit!" > GENES
    else
        echo "$hits" > GENES
    fi


>>>

output {
File blast_results = "~{samplename}.tblastn.tsv"
String blast_genes = read_string("GENES")
String blast_version = read_string("VERSION")
}

runtime {
    docker: "~{docker}"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible:  0
}
}


