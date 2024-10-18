version 1.0

task run_MASH {
    input {
        String sample_name
        File assembly
    }
    runtime{
        docker: 'staphb/mash:2.3'
    }

    command <<<
        date | tee DATE
        mash --version | tee VERSION

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
        String mash_version = read_string("VERSION")
        String mash_date = read_string("DATE")
        File mash_output = "~{sample_name}/mash_dist_output_sorted.txt"
    }

}