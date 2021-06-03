#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input          = "/home/AD/rbhuller/CODE/core/workflows/CNV/module_cnv/QC/samplesheet.csv" 
params.fasta          = "/home/AD/rbhuller/DATA/test_data/WGS_cell_line/genome/GRCh38_subseq_chr21.fasta"
params.annotationfile = "/home/AD/rbhuller/DATA/test_data/WGS_cell_line/illumina/txt/refFlat_hg38_15_without_chr.txt"


if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { INPUT_CHECK              } from './subworkflows/local/input_check'                   addParams( options: [:] )
include { FASTQC                   } from './modules/nf-core/software/fastqc/main'             addParams( options: [:] )
include { CUTADAPT                 } from './modules/nf-core/software/cutadapt/main'           addParams( options: [ 'args': '-q 20' ] )
include { FASTQC as CUTADAPTFASTQC } from './modules/nf-core/software/fastqc/main'             addParams( options: [:] )
include { BWAMEM2_INDEX            } from './modules/nf-core/software/bwamem2/index/main'      addParams( options: [:] )
include { BWAMEM2_MEM              } from './modules/nf-core/software/bwamem2/mem/main'        addParams( options: [:] )
include { BAM_STATS_SAMTOOLS       } from './subworkflows/nf-core/bam_stats_samtools'          addParams( options: [:] )
include { BAM_SORT_SAMTOOLS        } from './subworkflows/nf-core/bam_sort_samtools'           addParams( options: [:] )
include { BEDTOOLS_GENOMECOV       } from './modules/nf-core/software/bedtools/genomecov/main' addParams( options: [:] )
include { MARK_DUPLICATES_PICARD   } from './subworkflows/nf-core/mark_duplicates_picard'      addParams( options: [:] )
include { MULTIQC                  } from './modules/nf-core/software/multiqc/main'            addParams( options: [:] )
include { CNVKIT                   } from './modules/nf-core/software/cnvkit/main'             addParams( options: [ 'args': '--method wgs --output-reference reference.cnn' ] )

workflow test_QC {

    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .set { ch_fastq }
    

    FASTQC (ch_fastq)
    CUTADAPT (ch_fastq)
    CUTADAPTFASTQC (CUTADAPT.out.reads)
    BWAMEM2_INDEX (params.fasta)
    BWAMEM2_MEM (CUTADAPT.out.reads, BWAMEM2_INDEX.out.index)

    ch_bam = BWAMEM2_MEM.out.bam
    BAM_SORT_SAMTOOLS (ch_bam)
    
    BEDTOOLS_GENOMECOV (ch_bam)    

    bam = BAM_SORT_SAMTOOLS.out.bam
    MARK_DUPLICATES_PICARD (bam)

    MULTIQC (FASTQC.out.zip.mix(CUTADAPT.out.log, BAM_SORT_SAMTOOLS.out.all_results, MARK_DUPLICATES_PICARD.out.metrics).collect { it[1] })    


    lst = MARK_DUPLICATES_PICARD.out.bam.flatten().toList()
    rmlst = lst.remove(2)
    bamlst = lst.minus(rmlst)
    //bamlst.view()
    
    CNVKIT (bamlst, params.fasta, params.annotationfile)
}
