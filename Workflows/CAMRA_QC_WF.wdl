version 1.0

workflow assembly_qc {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run QC on a SINGLE pared-end illumina Assembly. Runs Mash, CheckM (inputs mash taxa), Merqury."
    }

    parameter_meta {
        sample_name     :   "Name of sample isolate"
        read1           :   "raw read 1 fastq.gz or fastq file"
        read2           :   "raw read 2 fastq.gz or fastq file"
        assembly        :   "assembly ins fasta or fasta.gz format"
        asm_size        :   "assembly size"
    }

    input{
        String sample_name
        File read1
        File read2
        File assembly
        String assembly_size
    }

    call run_MASH {
        input:
            sample_name = sample_name, 
            assembly = assembly
    }

    call run_entrez_direct {
        input:
            mashoutput = run_MASH.mash_output
    }

    call run_checkM {
        input:
            sample_name = sample_name,
            assembly = assembly,
            mash_genus = run_entrez_direct.mash_genus
    }

    call run_merqury {
        input:
            assembly = assembly,
            sample_name = sample_name,
            asm_size = assembly_size,
            read1 = read1,
            read2 = read2
    }

    output{
        String mash_ani = run_entrez_direct.mash_ani
        String mash_genus = run_entrez_direct.mash_genus
        String mash_species = run_entrez_direct.mash_species
        String mash_subspecies = run_entrez_direct.mash_subspecies
        String mash_taxaid = run_entrez_direct.mash_taxaid
        File checkm_output = run_checkM.checkm_output
        String checkm_markerlineage = run_checkM.checkm_markerlineage
        String checkm_completeness = run_checkM.checkm_completeness
        String checkm_contamination = run_checkM.checkm_contamination
        String checkm_heterogeneity = run_checkM.checkm_heterogeneity
        String merqury_qv = run_merqury.merqury_qv
        String merqury_comp = run_merqury.merqury_comp
        File merqury_qv_file = run_merqury.merqury_qv_file
        File merqury_completeness_file = run_merqury.merqury_completeness_file
    }


}

#Task runs MASH
task run_MASH {
    input {
        String sample_name
        File assembly
    }
    runtime{
        docker: 'staphb/mash:2.3'
    }

    command <<<
        #make sample dir and enter it
        mkdir ~{sample_name} && cd ~{sample_name}
        #Make a mash sketch of the assembly so that the mash can run faster
        mash sketch ~{assembly} -o ~{sample_name}
        # find the distance of the assembly to all the genomes in the reference database 
        mash dist /db/RefSeqSketchesDefaults.msh ~{sample_name}.msh > mash_dist_output.txt
        # sort the distance file by the mash value, we are intrested in the one with the smallest mash+
        sort -gk3 mash_dist_output.txt > mash_dist_output_sorted.txt
        # the first line of the sorted document will be our taxonomy
        # the line is compiosed of the following: tananomic id, bioprohject, biosample ID, ncbi refseq assembly, taxon (always genus, sometimes species, sometimes subspecies)
    >>>
    output{
        File mash_output = "~{sample_name}/mash_dist_output_sorted.txt"
    }

}

task run_entrez_direct{
    input {
        File mashoutput 
    }
    runtime {
        docker: 'danylmb/entrez-direct:latest'
    }
    command <<<
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




task run_checkM {
    input {
        String sample_name
        String mash_genus
        File assembly
    }
    runtime{
        docker: 'danylmb/checkm:latest'
    }
    command <<<
        export TMPDIR=/tmp
        mkdir assembly_dir
        
        if [[ "~{assembly}" == *.fasta.gz || "~{assembly}" == *.fa.gz ]]; then 
            gunzip -c ~{assembly} > assembly_dir/~{sample_name}.fasta 
        elif [[ "~{assembly}" == *.fasta || "~{assembly}" == *.fa ]]; then
            mv ~{assembly} assembly_dir/~{sample_name}.fasta
        fi

        export TMPDIR=/tmp
        if [[ -n "~{mash_genus}" ]]; then
            checkm taxonomy_wf genus ~{mash_genus} -t 4 -x fasta assembly_dir ~{sample_name} > checkm_quality_assessment.txt
        else
            checkm taxonomy_wf domain "Bacteria" -t 4 -x fasta assembly_dir  ~{sample_name} > checkm_quality_assessment.txt
        fi
        checkm_line=$(tail -n 3 "checkm_quality_assessment.txt" | head -n 1)
        read -r cm_ID cm_MarkerLineage cm1 cm2 cm3 cm4 cm5 cm6 cm7 cm8 cm9 cm10 cm_Completeness cm_Contamination cm_Heterogeneity <<< $checkm_line
        echo "$cm_MarkerLineage"
        echo "$cm_Completeness"
        echo "$cm_Contamination"
        echo "$cm_Heterogeneity"
    >>>
    output {
        File checkm_output = "checkm_quality_assessment.txt"
        Array[String] stdout_values = read_lines(stdout()) 
        String checkm_markerlineage = stdout_values[0]
        String checkm_completeness = stdout_values[1]
        String checkm_contamination = stdout_values[2]
        String checkm_heterogeneity = stdout_values[3]
    }
}



task run_merqury {
    input {
        File assembly
        File read1
        File read2
        String sample_name
        String asm_size
    }
    runtime{
        docker: 'miramastoras/merqury:latest'
        memory: "4G" #increasing memmory worked
    }
    command <<<
        # finding genome size and best kmer size

        total_length=~{asm_size} 
        best_k=$(best_k.sh $total_length) 
        best_k=$(echo "$best_k" | tail -n 1) 
        #Preparing meryl dbs
        meryl k=$best_k count output read1.meryl ~{read1}
        meryl k=$best_k count output read2.meryl ~{read2}

        meryl union-sum output ~{sample_name}.meryl read1.meryl read2.meryl 

        #Using Merqury
        mkdir merqury_output && cd merqury_output
        merqury.sh ../~{sample_name}.meryl ~{assembly} ~{sample_name} && echo "MERQURY DONE"
        #Outputs of Intrest
        qv_line=$(tail -n 3 ~{sample_name}.qv | head -n 1)
        qv_number=$(echo "$qv_line" | awk '{print $4}')
        echo $qv_number
        comp_line=$(tail -n 6 ~{sample_name}.completeness.stats | head -n 1)
        comp_number=$(echo "$comp_line" | awk '{print $5}')
        echo $comp_number

        >>>
    output {
        Array[String] stdout_values = read_lines(stdout()) 
        String merqury_qv = stdout_values[11]
        String merqury_comp = stdout_values[12]
        File merqury_qv_file = "merqury_output/~{sample_name}.qv"
        File merqury_completeness_file = "merqury_output/~{sample_name}.completeness.stats"
    }
}

