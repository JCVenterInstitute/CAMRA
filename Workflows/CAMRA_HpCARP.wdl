version 1.0

import "../Tasks/task_nanoplot.wdl" as nanoplot 
import "../Tasks/task_porechop.wdl" as porechop 
import "../Tasks/task_nanofilt.wdl" as nanofilt
import "../Tasks/task_minimap2.wdl" as minimap2
import "../Tasks/task_racon.wdl" as racon



workflow HpCARP {
    meta {
        author: "Daniella Matute"
        email: "dmatute@jcvi.org"
        description: "This WDL pipeline processes nanopore reads through quality checks, trimming, filtering, alignment, and consensus generation, outputting key metrics, processed reads, and consensus sequences."
    }
    
    parameter_meta {
        sample_name     :   "Name of sample isolate"
        ont_reads       :   "Raw ont reads, pooled reads of a single sample"
        reference_fasta :   "Reference SNPs"
        ref_seq         :   "Curated, non-redundant collection of reference sequences, representative of the central dogma, for each major organism."
    }

    input{
        String sample_name
        File ont_reads
        File reference_fasta
        File ref_seq
    }

    # nanoplot for basic QC metrics
    call nanoplot.run_nanoplot as nanoplot_raw {
        input:
            ont_reads = ont_reads,
            sample_name = sample_name
        }

    call porechop.run_porechop {
        input:
            ont_reads = ont_reads,
            sample_name = sample_name
        }

    call nanoplot.run_nanoplot as nanoplot_clean {
        input:
            ont_reads = run_porechop.porechop_reads,
            sample_name = sample_name
        }
    
    call nanofilt.run_nanofilt {
        input:
            porechop_reads = run_porechop.porechop_reads
        }

    call minimap2.run_minimap2 {
        input:
            fastq_file = run_nanofilt.nanofilt_fastq_output,
            reference_fasta = ref_seq
        }

    call racon.run_racon {
        input:
            minimap2_sorted_sam = run_minimap2.minimap2_sorted_output,
            nanofilt_fastq = run_nanofilt.nanofilt_fastq_output,
            reference_fasta = reference_fasta
        }

    output{
        # NANOPLOT_RAW
            File nanoplot_raw_html                  = nanoplot_raw.nanoplot_html
            File nanoplot_raw_tsv                   = nanoplot_raw.nanoplot_tsv
            String nanoplot_raw_num_reads              = nanoplot_raw.nanoplot_num_reads
            String nanoplot_raw_median_readlength    = nanoplot_raw.nanoplot_median_readlength
            String nanoplot_raw_mean_readlength      = nanoplot_raw.nanoplot_mean_readlength
            # Float nanoplot_raw_stdev_readlength     = nanoplot_raw.nanoplot_stdev_readlength
            String nanoplot_raw_n50                  = nanoplot_raw.nanoplot_n50
            String nanoplot_raw_mean_quality         = nanoplot_raw.nanoplot_mean_quality
            String nanoplot_raw_median_quality       = nanoplot_raw.nanoplot_median_quality
            String nanoplot_raw_version             = nanoplot_raw.nanoplot_version
            String nanoplot_raw_date                = nanoplot_raw.nanoplot_date

        # PORECHOP
            File porechop_reads                     = run_porechop.porechop_reads
            String porechop_version                 = run_porechop.porechop_version
            String porechop_date                    = run_porechop.porechop_date

        # NANOPLOT_CLEAN
            File nanoplot_clean_html                = nanoplot_clean.nanoplot_html
            File nanoplot_clean_tsv                 = nanoplot_clean.nanoplot_tsv
            String nanoplot_clean_num_reads            = nanoplot_clean.nanoplot_num_reads
            String nanoplot_clean_median_readlength  = nanoplot_clean.nanoplot_median_readlength
            String nanoplot_clean_mean_readlength    = nanoplot_clean.nanoplot_mean_readlength
            # Float nanoplot_clean_stdev_readlength   = nanoplot_clean.nanoplot_stdev_readlength
            String nanoplot_clean_n50                = nanoplot_clean.nanoplot_n50
            String nanoplot_clean_mean_quality       = nanoplot_clean.nanoplot_mean_quality
            String nanoplot_clean_median_quality     = nanoplot_clean.nanoplot_median_quality
            String nanoplot_clean_version           = nanoplot_clean.nanoplot_version
            String nanoplot_clean_date              = nanoplot_clean.nanoplot_date

        # NANOFILT
            File nanofilt_fastq_output              = run_nanofilt.nanofilt_fastq_output
            String nanofilt_version                 = run_nanofilt.nanofilt_version
        # MINIMAP2
            File minimap2_sorted_output             = run_minimap2.minimap2_sorted_output
            String minimap2_version                 = run_minimap2.minimap2_version
        # RACON
            File racon_consensus_output             = run_racon.racon_consensus_output
            String racon_version                    = run_racon.racon_version
    }
}
