version 1.0

task run_taxajoin {
    input {
        String genus
        String species
    }
    runtime{
        docker: "debian:stable-slim"
    }
    command <<<
        species_length=$(echo -n "~{species}" | wc -c)
        if [ "$species_length" -gt 2 ]; then
            # Concatenate "genus" and "species"
            organism="~{genus} ~{species}"
            echo "$organism"
        else
            organism="~{genus}"
            echo "Species probably invalid: ~{species}. Organism is set as: $organism"
        fi
    >>>

    output {
        
        String organism = read_string(stdout())
    }
}
