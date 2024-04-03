version 1.0

task run_entrez_direct{
    input {
        File mashoutput 
         #TODO add option to change docker version
    }
    runtime {
        docker: 'danylmb/entrez-direct:latest'
    }
    command <<<
        #TODO add versioning
        line=$(head -n 1 ~{mashoutput})
        taxa_id=$(echo "$line" | cut -d'-' -f3)
        scientific_name=$(efetch -db taxonomy -id "$taxa_id" -format xml | xtract -pattern Taxon -element ScientificName)
        read -r mash_genus mash_species mash_sub <<< "$scientific_name"
        mash_dist=$(echo "$line" | awk '{print $3}')
        ANI=$(echo "100 * ( 1 - $mash_dist )" | bc )
        echo "$ANI"

        empty=""
        # Define the thresholds
        threshold_1=97.5
        threshold_2=94
        threshold_3=80

        # Check conditions based on $ANI
        if (( $(bc -l <<< "$ANI >= $threshold_1") )); then
            :
        elif (( $(bc -l <<< "$ANI < $threshold_1 && $ANI >= $threshold_2") )); then
            mash_sub="$empty"
        elif (( $(bc -l <<< "$ANI < $threshold_2 && $ANI >= $threshold_3") )); then
            mash_species="$empty"
            mash_sub="$empty"
        else
            mash_genus="$empty"
            mash_species="$empty"
            mash_sub="$empty"
        fi

        # Echo variables
        echo "$mash_genus"
        echo "$mash_species"
        echo "$mash_sub"
        echo "$taxa_id"
    >>>
    output{
        Array[String] stdout_values = read_lines(stdout()) 
        String mash_ani = stdout_values[0]
        String mash_genus = stdout_values[1]
        String mash_species = stdout_values[2]
        String mash_subspecies = stdout_values[3]
        String mash_taxaid = stdout_values[4]
    }
}