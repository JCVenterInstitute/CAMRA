version 1.0

# TODO do you want to add the following tools? 
# plasmid finder or rgi or ariba or resfinder or virulence finder or abricate, fargene, deep arg
#ncbi-amrfinderplus, abricate, staramr, 
task run_AMRfinderPlus {

    meta {
    description: "AMRFinderPlus - Identify AMR genes and point mutations, and virulence and stress resistance genes in assembled bacterial nucleotide and protein sequence."
    gitrepository: "https://github.com/ncbi/amr"
    docker:"https://hub.docker.com/r/staphb/ncbi-amrfinderplus"
    cite: "Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984."
    } 
    input {
        File assembly
        String sample_name
        String organism
        String docker

    }
    runtime{
        docker: docker
    }

    command <<<

        date | tee date
        amrfinder --version |tee AMRFINDER_VERSION
        amrfinder --database_version 2>/dev/null | grep "Database version" | sed 's|Database version: ||' | tee AMRFINDER_DB_VERSION
        
        #TODO fix organism mapping
        # these are the organisms avaliable for Database version: 2024-01-31.1 when you run 'amrfinder -l' these might need to be manually updated is the data base version is increased
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
                echo "Either Gambit predicted taxon is not supported by NCBI-AMRFinderPlus or the user did not supply an organism as input."
                echo "Skipping the use of amrfinder --organism optional parameter.";;
        esac

        echo "amrfinder_organism is set to:" ${amrfinder_organism}
        #TODO Add plus flag
        # if amrfinder_organism variable is set, use --organism flag, otherwise do not use --organism flag
        if [[ -v amrfinder_organism ]] ; then
        # always use --plus flag, others may be left out if param is optional and not supplied 
            amrfinder --plus \
                --organism ${amrfinder_organism} \
                ~{'--name ' + sample_name} \
                ~{'--nucleotide ' + assembly} \
                ~{'-o ' + sample_name + '_amrfinder_all.tsv'} 
        else 
            echo "Either the organism (~{organism}) is not recognized by NCBI-AMRFinderPlus or the user did not supply an organism as input."
            echo "Skipping the use of amrfinder --organism optional parameter."
            # always use --plus flag, others may be left out if param is optional and not supplied 
            amrfinder --plus \
                ~{'--name ' + sample_name} \
                ~{'--nucleotide ' + assembly} \
                ~{'-o ' + sample_name + '_amrfinder_all.tsv'} 
            fi
    >>>
    output{
    }
}
