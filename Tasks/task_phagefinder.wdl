version 1.0

workflow annotation_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
    }

    parameter_meta {
        description: "Analysis of genome, AMR focused."
    }
    input {
        File annotation_txt
        File annotation_embl
        File assembly
        File annotation_protein
    }

    call run_PhageFinder {
        input:
            annotation_txt = annotation_txt,
            annotation_embl = annotation_embl,
            assembly = assembly,
            annotation_protein = annotation_protein
    }    

    output {
        String phagefinder_date = run_PhageFinder.phagefinder_date
    }

}


task run_PhageFinder {
    meta {
        description: "heuristic computer program written in PERL to identify prophage regions within bacterial genomes"
        gitrepository: "https://github.com/DanyMatute/Phage_Finder"
        publication: "https://pubmed.ncbi.nlm.nih.gov/17062630/"
    }

    input {
        File annotation_txt
        File annotation_embl
        File assembly
        File annotation_protein
    }

    runtime {
        docker: 'danylmb/phagefinder:22-11-2024'
    }

    command <<<
        date | tee DATE

        cp ~{assembly} ./sample.fna
        cp ~{annotation_protein} ./sample.faa

        python3 /opt/phage_finder_info-generator-BVBRC_annot.py ~{annotation_embl} ~{annotation_txt} . 

        /opt/Phage_Finder/bin/phage_finder.sh -s diamond ./sample
    >>>

    output {
        String phagefinder_date = read_string("DATE")
        #TODO finish outputs.

    }
}