version 1.0

import "../Tasks/task_checkm.wdl" as checkm 
import "../Tasks/task_mash.wdl" as mash 
import "../Tasks/task_entrezdirect.wdl" as entrezdirect
import "../Tasks/task_merqury.wdl" as merqury 
import "../Tasks/MLST.wdl" as mlst
import "../Tasks/task_quast.wdl" as quast
import "../Tasks/task_fastQC.wdl" as fastQC
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
        # TODO assembly size should be optional, if it is not provided then it should be calculated another way, 
        # like by a quast or busco analysis, does checkm produce it??
        String assembly_size 
        Int? min_contigs
    }

    # call mash.run_MASH {
    #     input:
    #         sample_name = sample_name, 
    #         assembly = assembly
    # }
    
    # call mlst.run_MLST {
    #     input:
    #         assembly = assembly,
    #         sample_name = sample_name
    # }

    # call entrezdirect.run_entrez_direct {
    #     input:
    #         mashoutput = run_MASH.mash_output
    # }

    # call checkm.run_checkM {
    #     input:
    #         sample_name = sample_name,
    #         assembly = assembly,
    #         mash_genus = run_entrez_direct.mash_genus
    # }

    # call merqury.run_merqury {
    #     input:
    #         assembly = assembly,
    #         sample_name = sample_name,
    #         asm_size = assembly_size,
    #         read1 = read1,
    #         read2 = read2
    # }

    # call quast.run_Quast {
    #     input:
    #         assembly = assembly,
    #         min_contigs = min_contigs
    # }

    call fastQC.run_fastQC {
        input:
            read1 = read1,
            read2 = read2,
    }

    output{
        File fastQC_R1_data = run_fastQC.fastQC_R1_data
        File fastQC_R2_data = run_fastQC.fastQC_R2_data
        File fastqc_R1_report = run_fastQC.fastQC_R1_html
        File fastqc_R2_report = run_fastQC.fastQC_R2_html
        File fastQC_R1_summary = run_fastQC.fastQC_R1_summary
        File fastQC_R2_summary = run_fastQC.fastQC_R2_summary

        String fastQC_R1_PassWarnFail = run_fastQC.fastQC_R1_PassWarnFail
        String fastQC_R2_PassWarnFail = run_fastQC.fastQC_R2_PassWarnFail

        String fastQC_R1_total_sequences = run_fastQC.fastQC_R1_total_sequences
        String fastQC_R1_total_bases = run_fastQC.fastQC_R1_total_bases
        String fastQC_R1_poor_quality = run_fastQC.fastQC_R1_poor_quality
        String fastQC_R1_sequence_length = run_fastQC.fastQC_R1_sequence_length
        String fastQC_R1_gc_content = run_fastQC.fastQC_R1_gc_content
        String fastQC_R2_total_sequences = run_fastQC.fastQC_R2_total_sequences
        String fastQC_R2_total_bases = run_fastQC.fastQC_R2_total_bases
        String fastQC_R2_poor_quality = run_fastQC.fastQC_R2_poor_quality
        String fastQC_R2_sequence_length = run_fastQC.fastQC_R2_sequence_length
        String fastQC_R2_gc_content = run_fastQC.fastQC_R2_gc_content
    
        # File quast_report = run_Quast.quast_report
        # Int quast_largest_contig_value = run_Quast.quast_contig_largest
        # Int quast_total_length_value = run_Quast.quast_total_length
        # Int quast_n50_value = run_Quast.quast_N50
        # Int quast_n90_value = run_Quast.quast_N90
        # Int quast_l50_value = run_Quast.quast_L50
        # Int quast_l90_value = run_Quast.quast_L90

        # String mash_ani = run_entrez_direct.mash_ani
        # String mash_genus = run_entrez_direct.mash_genus
        # String mash_species = run_entrez_direct.mash_species
        # String mash_subspecies = run_entrez_direct.mash_subspecies
        # String mash_taxaid = run_entrez_direct.mash_taxaid

        # File checkm_output = run_checkM.checkm_output
        # String checkm_markerlineage = run_checkM.checkm_markerlineage
        # String checkm_completeness = run_checkM.checkm_completeness
        # String checkm_contamination = run_checkM.checkm_contamination
        # String checkm_heterogeneity = run_checkM.checkm_heterogeneity

        # String merqury_qv = run_merqury.merqury_qv
        # String merqury_comp = run_merqury.merqury_comp
        # File merqury_qv_file = run_merqury.merqury_qv_file
        # File merqury_completeness_file = run_merqury.merqury_completeness_file

        # String mlst_scheme = run_MLST.tsMLST_scheme 
        # String mlst_seqtype = run_MLST.tsMLST_seqtype 
        # String mlst_alleles = run_MLST.tsMLST_alleles

    }


}
