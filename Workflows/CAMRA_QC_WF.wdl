version 1.0

import "../Tasks/task_checkm.wdl" as checkm 
import "../Tasks/task_mash.wdl" as mash 
import "../Tasks/task_entrezdirect.wdl" as entrezdirect
import "../Tasks/task_merqury.wdl" as merqury 
import "../Tasks/task_mlst.wdl" as mlst
import "../Tasks/task_quast.wdl" as quast
import "../Tasks/task_fastQC.wdl" as fastQC

workflow assembly_qc {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "Run QC on a SINGLE pared-end illumina Assembly. Runs Mash, CheckM (inputs mash taxa), Merqury."
    }

    parameter_meta {
        sample_name     :   "Name of sample isolate"
        read1           :   "Raw read 1 fastq.gz or fastq file"
        read2           :   "Raw read 2 fastq.gz or fastq file"
        assembly        :   "Assembly in fasta or fasta.gz format"
        pubmlst_DB      :   "Optional pubmlst_DB for mlst mapping"
    }

    input{
        String sample_name
        File read1
        File read2
        File assembly
        File pubmlst_DB = "None"
    }

    call mash.run_MASH {
        input:
            sample_name = sample_name, 
            assembly = assembly
    }
    
    if ( "~{pubmlst_DB}" == "None") {
        call mlst.run_MLST {
            input:
                assembly = assembly
        }
    }
    
    if ( "~{pubmlst_DB}" != "None" ) {
        call mlst.run_MLST_pubmlst_DB {
            input: 
                assembly = assembly,
                pubmlst_DB = pubmlst_DB
        }
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

    call quast.run_Quast {
        input:
            assembly = assembly
            #min_contigs = min_contigs
    }

    call merqury.run_merqury {
        input:
            assembly = assembly,
            asm_size = run_Quast.quast_total_length,
            read1 = read1,
            read2 = read2
    }

    call fastQC.run_fastQC {
        input:
            read1 = read1,
            read2 = read2
    }

    output{
        # Quast
        String quast_version = run_Quast.quast_version
        String quast_date = run_Quast.quast_date
        File quast_report = run_Quast.quast_report
        Int quast_largest_contig = run_Quast.quast_contig_largest
        Int quast_asm_length = run_Quast.quast_total_length
        Int quast_N50 = run_Quast.quast_N50
        Int quast_N90 = run_Quast.quast_N90
        Int quast_L50 = run_Quast.quast_L50
        Int quast_L90 = run_Quast.quast_L90

        #MASH
        String mash_version = run_MASH.mash_version
        String mash_date = run_MASH.mash_date
        String mash_ani = run_entrez_direct.mash_ani
        String mash_genus = run_entrez_direct.mash_genus
        String mash_species = run_entrez_direct.mash_species
        String mash_subspecies = run_entrez_direct.mash_subspecies
        String mash_taxaid = run_entrez_direct.mash_taxaid

        #CheckM
        String checkm_version = run_checkM.checkm_version
        String checkm_date = run_checkM.checkm_date
        File checkm_output = run_checkM.checkm_output
        String checkm_markerlineage = run_checkM.checkm_markerlineage
        String checkm_completeness = run_checkM.checkm_completeness
        String checkm_contamination = run_checkM.checkm_contamination
        String checkm_heterogeneity = run_checkM.checkm_heterogeneity

        #Merqury
        String merqury_qv = run_merqury.merqury_qv
        String merqury_completeness = run_merqury.merqury_comp
        File merqury_qv_file = run_merqury.merqury_qv_file
        File merqury_completeness_file = run_merqury.merqury_completeness_file
        String merqury_version = run_merqury.merqury_version
        String merqury_date = run_merqury.merqury_date

        #MLST
        String? tsMLST_version = run_MLST.tsMLST_version
        String? tsMLST_date = run_MLST.tsMLST_date
        File? tsMLST_tsv_output = run_MLST.tsMLST_tsv_output
        String? tsMLST_scheme = run_MLST.tsMLST_scheme 
        String? tsMLST_seqtype = run_MLST.tsMLST_seqtype 
        String? tsMLST_alleles = run_MLST.tsMLST_alleles

        #MLST + pubmlst_DB
        String? tsMLSTdb_version = run_MLST_pubmlst_DB.tsMLST_version
        String? tsMLSTdb_date = run_MLST_pubmlst_DB.tsMLST_date
        File? tsMLSTdb_tsv_output = run_MLST_pubmlst_DB.tsMLST_tsv_output
        String? tsMLSTdb_scheme = run_MLST_pubmlst_DB.tsMLST_scheme 
        String? tsMLSTdb_seqtype = run_MLST_pubmlst_DB.tsMLST_seqtype 
        String? tsMLSTdb_alleles = run_MLST_pubmlst_DB.tsMLST_alleles

        #FastQC
        File fastQC_R1_html = run_fastQC.fastQC_R1_html
        File fastQC_R2_html = run_fastQC.fastQC_R2_html

        String fastQC_R1_PassWarnFail = run_fastQC.fastQC_R1_PassWarnFail
        String fastQC_R2_PassWarnFail = run_fastQC.fastQC_R2_PassWarnFail

        String fastQC_R1_total_reads = run_fastQC.fastQC_R1_total_reads
        String fastQC_R1_total_bases = run_fastQC.fastQC_R1_total_bases
        String fastQC_R1_read_length = run_fastQC.fastQC_R1_read_length 
        String fastQC_R1_gc_content = run_fastQC.fastQC_R1_gc_content

        String fastQC_R2_total_reads = run_fastQC.fastQC_R2_total_reads
        String fastQC_R2_total_bases = run_fastQC.fastQC_R2_total_bases
        String fastQC_R2_read_length = run_fastQC.fastQC_R2_read_length
        String fastQC_R2_gc_content = run_fastQC.fastQC_R2_gc_content
        
        String fastQC_version = run_fastQC.fastQC_version
        String fastQC_date = run_fastQC.fastQC_date
    }


}
