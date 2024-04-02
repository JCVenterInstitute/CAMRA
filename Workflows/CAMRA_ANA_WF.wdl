version 1.0
workflow assembly_qc {
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
        String? dockerversion_AMRfinder
    }
}
# plasmid finder or rgi or ariba or resfinder or virulence finder or abricate, fargene, deep arg
#ncbi-amrfinderplus, abricate, staramr, 
task run_AMRfinderPlus {
    input {
        File assembly
        String sample_name
        String organism
        String docker = "staphb/ncbi-amrfinderplus:3.12.8"

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


    '''Available --organism options: Acinetobacter_baumannii, Burkholderia_cepacia, Burkholderia_pseudomallei, Campylobacter, 
    Citrobacter_freundii, Clostridioides_difficile, Enterobacter_asburiae, Enterobacter_cloacae, Enterococcus_faecalis, 
    Enterococcus_faecium, Escherichia, Klebsiella_oxytoca, Klebsiella_pneumoniae, Neisseria_gonorrhoeae, Neisseria_meningitidis, 
    Pseudomonas_aeruginosa, Salmonella, Serratia_marcescens, Staphylococcus_aureus, Staphylococcus_pseudintermedius, 
    Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae, Vibrio_parahaemolyticus, Vibrio_vulnificus'''
        # if amrfinder_organism variable is set, use --organism flag, otherwise do not use --organism flag
    if [[ -v amrfinder_organism ]] ; then
      # always use --plus flag, others may be left out if param is optional and not supplied 
      amrfinder --plus \
        --organism ${amrfinder_organism} \
        ~{'--name ' + samplename} \
        ~{'--nucleotide ' + assembly} \
        ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
        ~{'--threads ' + cpu} \
        ~{'--coverage_min ' + mincov} \
        ~{'--ident_min ' + minid}
    else 
      echo "Either the organism (~{organism}) is not recognized by NCBI-AMRFinderPlus or the user did not supply an organism as input."
      echo "Skipping the use of amrfinder --organism optional parameter."
      # always use --plus flag, others may be left out if param is optional and not supplied 
      amrfinder --plus \
        ~{'--name ' + samplename} \
        ~{'--nucleotide ' + assembly} \
        ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
        ~{'--threads ' + cpu} \
        ~{'--coverage_min ' + mincov} \
        ~{'--ident_min ' + minid}
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