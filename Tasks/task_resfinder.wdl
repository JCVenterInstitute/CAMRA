version 1.0

task run_ResFinder {
    meta {
    description: "ResFinder identifies (BLAST & KMA) acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria."
    repository:"https://bitbucket.org/genomicepidemiology/resfinder/src/4.5.0/"
    cite: "ResFinder 4.0 for predictions of phenotypes from genotypes. Bortolaia V, Kaas RS, Ruppe E, Roberts MC, Schwarz S, Cattoir V, Philippon A, Allesoe RL, Rebelo AR, Florensa AR, Fagelhauer L, Chakraborty T, Neumann B, Werner G, Bender JK, Stingl K, Nguyen M, Coppens J, Xavier BB, Malhotra-Kumar S, Westh H, Pinholt M, Anjum MF, Duggett NA, Kempf I, Nykasenoja S, Olkkola S, Wieczorek K, Amaro A, Clemente L, Mossong J, Losch S, Ragimbeau C, Lund O, Aarestrup FM. Journal of Antimicrobial Chemotherapy. 2020 Aug 11. PMID: 32780112 doi: 10.1093/jac/dkaa345"
    } 

    runtime{
        docker: 'danylmb/resfinder:v4.5.0'
        memory: "8 GB"
        cpu: 4
        disks: "local-disk " + 100 + " SSD"
        disk: 100 + " GB"
    }

    input {
        File assembly
        String? organism
        Float? min_cov = 0.6 # Min coverage 
        Float? threshold = 0.9 # Min identity
    }

    command <<<

        python3 -m resfinder --version | tee RESFINDER_VERSION
        # Versions obtianed from the dockerfile
        echo "1.4.14" | tee KMA_VERSION 
        echo "resfinder_db=2.3.1, pointfinder_db=4.1.0, disinfinder_db=2.0.1" | tee DB_VERSION
        date | tee DATE
        

        unzip_and_move() {
            local file="$1"
            local output="$2"
            if [[ "$file" == *.gz ]]; then
                echo "We need to unzip $file."
                gunzip -c "$file" > "$output" && echo "    > unzip & mv completed"
            else
                echo "We do not need to unzip $file."
                mv "$file" "$output" && echo "    > mv successful"
            fi
        }

        unzip_and_move "~{assembly}" "./assembly.fasta"

        # Make directories that will store the resfinder run wiht the assembly file and the read file 
        mkdir assembly_output

        case "~{organism}" in
            *"Acinetobacter"*"baumannii"*)
                species="Acinetobacter_baumannii";;
            *"Burkholderia"*"cepacia"*)
                species="Burkholderia_cepacia";;
            *"Burkholderia"*"pseudomallei"*)
                species="Burkholderia_pseudomallei";;
            *"Campylobacter"*"coli"* | *"Campylobacter"*"jejuni"*)
                species="Campylobacter";;
            *"Citrobacter"*"freundii"*)
                species="Citrobacter_freundii";;
            *"Clostridioides"*"difficile"*)
                species="Clostridioides_difficile";;
            *"Enterobacter"*"asburiae"*)
                species="Enterobacter_asburiae";;
            *"Enterobacter"*"asburiae"*)
                species="Enterobacter_asburiae";;
            *"Enterococcus"*"faecalis"*)
                species="Enterococcus_faecalis";;
            *"Enterococcus"*"faecium"* | *"Enterococcus"*"hirae"*)
                species="Enterococcus_faecium";;
            *"Escherichia"* | *"Shigella"*)
                species="Escherichia";;
            *"Klebsiella"*"oxytoca"*)
                species="Klebsiella_oxytoca";;
            *"Klebsiella"*"pneumoniae"*)
                species="Klebsiella_pneumoniae";;
            *"Neisseria"*"gonorrhoeae"*)
                species="Neisseria_gonorrhoeae";;
            *"Neisseria"*"meningitidis"*)
                species="Neisseria_meningitidis";;
            *"Pseudomonas"*"aeruginosa"*)
                species="Pseudomonas_aeruginosa";;
            *"Salmonella"*)
                species="Salmonella";;
            *"Serratia"*"marcescens"*)
                species="Serratia_marcescens";;
            *"Staphylococcus"*"aureus"*)
                species="Staphylococcus_aureus";;
            *"Staphylococcus"*"pseudintermedius"*)
                species="Staphylococcus_pseudintermedius";;
            *"Streptococcus"*"agalactiae"*)
                species="Streptococcus_agalactiae";;
            *"Streptococcus"*"pneumoniae"* | *"Streptococcus"*"mitis"*)
                species="Streptococcus_pneumoniae";;
            *"Streptococcus"*"pyogenes"*)
                species="Streptococcus_pyogenes";;
            *"Vibrio"*"cholerae"*)
                species="Vibrio_cholerae";;
            *"Vibrio"*"parahaemolyticus"*)
                species="Vibrio_parahaemolyticus";;
            *"Vibrio"*"vulnificus"*)
                species="Vibrio_vulnificus";;
            *)
                echo "species is not mapped.";;
        esac

        echo "species is set to:" ${species}

        # Run the assembly file on resfinder, only does BLAST
        python3 -m resfinder --inputfasta assembly.fasta --species species  --disinfectant --acquired --point --ignore_missing_species --outputPath assembly_output  ~{'--min_cov ' + min_cov}  ~{'--threshold ' + threshold} 


        rm assembly.fasta

        mv assembly_output/ResFinder_results_tab.txt assembly_output/resfinder_asm_results_tab.txt

        >>>
    output {
        String resfinder_version = read_string("RESFINDER_VERSION")
        String resfinder_kma_version = read_string("KMA_VERSION")
        String resfinder_db_version = read_string("DB_VERSION")
        String resfinder_date = read_string("DATE")

        File resfider_asm_output = "assembly_output/resfinder_asm_results_tab.txt"
        # File? resfinder_read_output = "read_output/resfinder_read_results_tab.txt"
        File resfinder_asm_hits = "assembly_output/ResFinder_Hit_in_genome_seq.fsa"
        # File? resfinder_read_hits = "read_output/ResFinder_Hit_in_genome_seq.fsa"
        File resfinder_asm_argseq = "assembly_output/ResFinder_Resistance_gene_seq.fsa"
        # File? resfinder_read_argseq = "read_output/ResFinder_Resistance_gene_seq.fsa"
    }
}