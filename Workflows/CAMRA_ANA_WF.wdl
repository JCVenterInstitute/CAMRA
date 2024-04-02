version 1.0
workflow assembly_analysis   {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run Analysis on a SINGLE Assembly."
    }

    parameter_meta {
        description: "Analysis of genomes"
    }
    input {
        String sample_name
        # File read1
        # File read2
        File assembly
        String organism 
        String docker
    }

    call run_AMRfinderPlus {
        input:
            assembly = assembly,
            sample_name = sample_name,
            organism = organism,
            docker = docker  


    }
}
# plasmid finder or rgi or ariba or resfinder or virulence finder or abricate, fargene, deep arg
#ncbi-amrfinderplus, abricate, staramr, 
task run_AMRfinderPlus {
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
        
        # these are the organisms avaliable for Database version: 2024-01-31.1 when you run 'amrfinder -l' these might need to be manually updated is the data base version is increased
        case "$organism" in
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

# #planton, 
# task run_PlasmidFinder {
#     input {

#     }
#     runtime{
#         docker: docker
#     }

#     command <<<
#     >>>
#     output{

#     }

# }

# #with Mash or Kraken?
# task run_Species_Typing {
#     input {

#     }
#     runtime{
#         docker: 
#     }

#     command <<<
#     >>>
#     output{

#     }

# }

# task run_Phage_Finder {
#     input {

#     }
#     runtime{
#         docker: 
#     }

#     command <<<
#     >>>
#     output{

#     }

# }