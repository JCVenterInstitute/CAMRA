version 1.0

task run_AMRfinderPlus {

    meta {
    description: "AMRFinderPlus - Identify AMR genes and point mutations, and virulence and stress resistance genes in assembled bacterial nucleotide and protein sequence."
    gitrepository: "https://github.com/ncbi/amr"
    docker:"https://hub.docker.com/r/staphb/ncbi-amrfinderplus"
    cite: "Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984."
    } 
    input {
        File assembly
        String organism
    }
    runtime{
        docker: "staphb/ncbi-amrfinderplus:3.12.8-2024-01-31.1"
    }

    command <<<

        date | tee DATE
        amrfinder --version |tee VERSION
        amrfinder --database_version 2>/dev/null | grep "Database version" | sed 's|Database version: ||' | tee DB_VERSION
        
        # When you run 'amrfinder -l' you get a list of the available organisms. Here, the sample's organism is aligned to its match within the database if available. 
            # The list was created with the databse wersion 2024-01-31.1 on April 2nd 2024. 
            # The list will be updated in subsequent versions of this wdl. Users are welcome to update it as well.  
        case "~{organism}" in
            *"Acinetobacter"*"baumannii"*)
                amrfinder_organism="Acinetobacter_baumannii";;
            *"Burkholderia"*"cepacia"*)
                amrfinder_organism="Burkholderia_cepacia";;
            *"Burkholderia"*"pseudomallei"*)
                amrfinder_organism="Burkholderia_pseudomallei";;
            *"Campylobacter"*"coli"* | *"Campylobacter"*"jejuni"*)
                amrfinder_organism="Campylobacter";;
            *"Citrobacter"*"freundii"*)
                amrfinder_organism="Citrobacter_freundii";;
            *"Clostridioides"*"difficile"*)
                amrfinder_organism="Clostridioides_difficile";;
            *"Enterobacter"*"asburiae"*)
                amrfinder_organism="Enterobacter_asburiae";;
            *"Enterobacter"*"asburiae"*)
                amrfinder_organism="Enterobacter_asburiae";;
            *"Enterococcus"*"faecalis"*)
                amrfinder_organism="Enterococcus_faecalis";;
            *"Enterococcus"*"faecium"* | *"Enterococcus"*"hirae"*)
                amrfinder_organism="Enterococcus_faecium";;
            *"Escherichia"* | *"Shigella"*)
                amrfinder_organism="Escherichia";;
            *"Klebsiella"*"oxytoca"*)
                amrfinder_organism="Klebsiella_oxytoca";;
            *"Klebsiella"*"pneumoniae"*)
                amrfinder_organism="Klebsiella_pneumoniae";;
            *"Neisseria"*"gonorrhoeae"*)
                amrfinder_organism="Neisseria_gonorrhoeae";;
            *"Neisseria"*"meningitidis"*)
                amrfinder_organism="Neisseria_meningitidis";;
            *"Pseudomonas"*"aeruginosa"*)
                amrfinder_organism="Pseudomonas_aeruginosa";;
            *"Salmonella"*)
                amrfinder_organism="Salmonella";;
            *"Serratia"*"marcescens"*)
                amrfinder_organism="Serratia_marcescens";;
            *"Staphylococcus"*"aureus"*)
                amrfinder_organism="Staphylococcus_aureus";;
            *"Staphylococcus"*"pseudintermedius"*)
                amrfinder_organism="Staphylococcus_pseudintermedius";;
            *"Streptococcus"*"agalactiae"*)
                amrfinder_organism="Streptococcus_agalactiae";;
            *"Streptococcus"*"pneumoniae"* | *"Streptococcus"*"mitis"*)
                amrfinder_organism="Streptococcus_pneumoniae";;
            *"Streptococcus"*"pyogenes"*)
                amrfinder_organism="Streptococcus_pyogenes";;
            *"Vibrio"*"cholerae"*)
                amrfinder_organism="Vibrio_cholerae";;
            *"Vibrio"*"parahaemolyticus"*)
                amrfinder_organism="Vibrio_parahaemolyticus";;
            *"Vibrio"*"vulnificus"*)
                amrfinder_organism="Vibrio_vulnificus";;
            *)
                echo "amrfinder_organism is not mapped.";;
        esac

        echo "amrfinder_organism is set to:" ${amrfinder_organism}


        # if amrfinder_organism variable is set, use --organism flag, otherwise do not use --organism flag
        if [[ -v amrfinder_organism ]] ; then
            echo "Running AMRFinder+ WITH amrfinder_organism."
            # always use --plus flag, others may be left out if param is optional and not supplied 
            amrfinder --plus \
                --organism ${amrfinder_organism} \
                ~{'--nucleotide ' + assembly} \
                >  amrfinderplus_output_all.txt

        else 
            echo "Running AMRFinder+ WITHOUT amrfinder_organism."
            # always use --plus flag, others may be left out if param is optional and not supplied 
            amrfinder --plus \
                ~{'--nucleotide ' + assembly} \
                >  amrfinderplus_output_all.txt
            fi

        # Element Type possibilities: AMR, STRESS, and VIRULENCE 
        # create headers for 3 output files; tee to 3 files and redirect STDOUT to dev null so it doesn't print to log file
        head -n 1 amrfinderplus_output_all.txt | tee amrfinderplus_stress.tsv amrfinderplus_virulence.tsv amrfinderplus_amr.tsv >/dev/null
        # looks for all rows with STRESS, AMR, or VIRULENCE and append to TSVs
        grep -F 'STRESS' amrfinderplus_output_all.txt >> amrfinderplus_stress.tsv || true
        grep -F 'VIRULENCE' amrfinderplus_output_all.txt >> amrfinderplus_virulence.tsv || true
        # || true is so that the final grep exits with code 0, preventing failures
        grep -F 'AMR' amrfinderplus_output_all.txt >> amrfinderplus_amr.tsv || true
    >>>
    output{

        String amrinder_version = read_string("VERSION")
        String amrfinder_db_version = read_string("DB_VERSION")
        String amrfinder_date = read_string("DATE")

        File amrfinder_all_output = "amrfinderplus_output_all.txt"
        File amrfinder_stress_output = "amrfinderplus_stress.tsv"
        File amrfinder_virulence_output = "amrfinderplus_virulence.tsv"
        File amrfinder_amr_output = "amrfinderplus_amr.tsv"
    }
}

task run_Query_Blastn {

    meta {
    description: "Blastn genome against query"
    gitrepository: "https://github.com/ncbi/amr"
    docker:"https://hub.docker.com/r/staphb/ncbi-amrfinderplus"
    cite: "Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984."
    } 
    input {
        File assembly
        File? query
    }
    
    runtime{
        docker: "staphb/ncbi-amrfinderplus:3.12.8-2024-01-31.1"
    }

    command <<<

        # Unzips the assembly file if nessesary
        if [[ "~{assembly}" == *.fasta.gz || "~{assembly}" == *.fa.gz || "~{assembly}" == *.fna.gz ]]; then 
            gunzip -c ~{assembly} > /data/assembly.fasta 
        elif [[ "~{assembly}" == *.fasta || "~{assembly}" == *.fa || "~{assembly}" == *.fna.gz ]]; then
            mv ~{assembly} /data/assembly.fasta
        fi

        # makeblastdb sets up a structured, searchable version of your genome sequence so that BLAST can efficiently find and match your target genes.
        makeblastdb -in /data/assembly.fasta -dbtype nucl -out genome_db

        # BLASTn command that detects variations and closely related sequences effectively
        blastn -query  ~{query} -db genome_db -out results.out -outfmt 6 -evalue 1e-10 -reward 2 -penalty -3 -word_size 7 -max_target_seqs 10 -dust no
        
        # sort results.out by the bitscore
        sort -k 12,12nr -t$'\t' results.out -o blast_results.txt
    >>>
    output{
        File blastn_output = "blast_results.txt" 
    }
}
