version 1.0

import "../Tasks/task_checkm.wdl" as checkm 
import "../Tasks/task_mash.wdl" as mash 
import "../Tasks/task_entrezdirect.wdl" as entrezdirect
import "../Tasks/task_merqury.wdl" as merqury 
import "../Tasks/MLST.wdl" as mlst
#TODO add quast and busco and or other qc tools

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
        #TODO assembly size should be optional, if it is not provided then it should be calculated another way
        String assembly_size 
    }

    call mash.run_MASH {
        input:
            sample_name = sample_name, 
            assembly = assembly
    }
    
    call mlst.run_MLST {
        input:
            assembly = assembly,
            sample_name = sample_name
    }

    call entrezdirect.run_entrez_direct {
        input:
            mashoutput = run_MASH.mash_output
    }

    call checkm.run_checkM {
        input:
            sample_name = sample_name,
            assembly = assembly,
            mash_genus = run_entrez_direct.mash_genus
    }

    call merqury.run_merqury {
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

        String mlst_scheme = run_MLST.tsMLST_scheme 
        String mlst_seqtype = run_MLST.tsMLST_seqtype 
        String mlst_alleles = run_MLST.tsMLST_alleles

    }


}
