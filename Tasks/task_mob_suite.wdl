version 1.0

task run_mob_suite {
    meta {
        description: "Find Plasmids using mob suite"
        version: "1.0"
        dockerhub: "bessonov/mob_suite:3.0.3"
    }

    input {
        File assembly,
        String SampleID
    }

    runtime {
        docker: "kbessonov/mob_suite:3.0.3"
    }

    command <<<
        if [[ "~{assembly}" == *.fasta.gz || "~{assembly}" == *.fa.gz ]]; then 
            gunzip -c ~{assembly} > ~{SampleID}.fasta 
        elif [[ "~{assembly}" == *.fasta || "~{assembly}" == *.fa ]]; then
            mv ~{assembly}  ~{SampleID}.fasta
        fi

        mob_recon -i ~{SampleID}.fasta -o ./mob/
        if [ ! -f "./mob/mobtyper_results.txt" ]; then
             touch ./mob/mobtyper_results.txt
        fi
        mob_typer -i assembly.fasta -o ./mob/Typing.out 
    >>>

    output {
        File mob_contig_out = "mob/contig_report.txt"
        File mob_plasmid_out = "mob/mobtyper_results.txt"
        File mob_typer_out = "mob/Typing.out"
    }
}
