version 1.0

task run_checkM {
    input {
        String sample_name
        String mash_genus
        File assembly
        #TODO add option to change docker version
    }
    runtime{
        docker: 'danylmb/checkm:latest'
    }
    command <<<
        export TMPDIR=/tmp
        mkdir assembly_dir
        #TODO add versioning
        
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