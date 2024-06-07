version 1.1

task RunRGI { 
    input {
        File output_sorted_length_100
        File card_json
        File kmer_db
        File amr_kmer_db
        File wildcard_data 
        File wildcard_index
        String docker_image_id
    }
    command <<<
        set -exuo pipefail
        time rgi load \
            -i "~{card_json}" \
            --wildcard_annotation "~{wildcard_data}" \
            --wildcard_version 4.0.0 \
            --wildcard_index "~{wildcard_index}" \
            --kmer_database "~{kmer_db}" \
            --amr_kmers "~{amr_kmer_db}" \
            --kmer_size 61
        rgi kmer_query --bwt -k 61 -i "~{output_sorted_length_100}" --output sr_species_report
    >>>
    output {
        Array[File]+ output_kmer_bwt = glob("sr_species_report*")
        File sr_species_allele = "sr_species_report_61mer_analysis.allele.txt"
        File kma_species_calling = "sr_species_report_61mer_analysis.gene.txt"
    }
    runtime {
        docker: docker_image_id
    }

}