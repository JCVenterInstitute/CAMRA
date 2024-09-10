version 1.0

task run_taxajoin {
    input {
        String genus
        String species
    }
    runtime{
        docker: "staphb/ncbi-amrfinderplus:3.12.8-2024-01-31.1"
    }
    command <<<
        species_length=$(echo -n "~{species}" | wc -c)
        if [ "$species_length" -gt 3 ]; then
            # Concatenate "genus" and "species"
            organism="~{genus} ~{species}"
            echo "Organism is set as: $organism"
        else
            organism="~{genus}"
            echo "Species probably invalid: ~{species}. Organism is set as: $organism"
        fi
    >>>

    output {
        
        String organism = read_string(stdout())
    }
}