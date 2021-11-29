#!/usr/bin/env nextflow
// This file defines individual processes (separated for portability)

// perform BAM merging with samtools
process "samtools_merge" {

    label "low"
    label "finish"

    input:
    path bams
    // eg. [sample1.bam, sample2.bam, ... sampleN.bam]

    output:
    path "merged/merged.bam"

    script:
    """
    mkdir merged
    samtools merge merged/merged.bam ${bams}
    """
}


// mark PCR duplicates with Picard MarkDuplicates
process "Picard_MarkDuplicates" {

    label "low"
    label "finish"

    publishDir "${params.output}", pattern: "bam/*.{bam,txt}", mode: 'copy'
    publishDir "${params.output}", pattern: "bam/logs/markDups.log", mode: 'move'

    input:
    path bam

    output:
    path "bam/*.bam"
    path "bam/*.txt"
    path "bam/logs/markDups.log"

    when:
    params.markDups

    script:
    """
    mkdir tmp bam bam/logs
    picard MarkDuplicates TMP_DIR=tmp \\
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=\$(ulimit -n) \\
    VALIDATION_STRINGENCY=LENIENT \\
    I=${bam} O=bam/markDups.bam \\
    M=bam/duplicates.txt > bam/logs/markDups.log 2>&1
    """
}



// perform QC of alignments
process "bamQC" {

    label "low"
    label "ignore"

    publishDir "${params.output}", pattern: "bam/*.{txt,pdf}", mode: 'move'
    publishDir "${params.output}", pattern: "bam/logs/bamQC.log", mode: 'move'

    input:
    path bam

    output:
    path "bam/*.{txt,pdf}"
    path "bam/logs/bamQC.log"

    when:
    params.bamQC

    script:
    """
    mkdir bam bam/logs
    qualimap bamqc -bam ${bam} -outdir bam -outformat pdf > bam/logs/bamQC.log 2>&1
    """
}



// joint variant calling with Freebayes
process "Freebayes" {

    label "low"
    label "finish"

    publishDir "${params.output}", pattern: "vcf/raw.vcf", mode: params.minQual > 0 ? 'copy' : 'move'

    input:
    path bam
    path fasta
    path fai

    output:
    path "vcf/raw.vcf"

    script:
    """
    mkdir vcf
    freebayes -f ${fasta} ${bam} > vcf/raw.vcf
    """
}



// perform quality filtering of VCF files with bcftools
process "bcftools" {

    label "low"
    label "finish"

    publishDir "${params.output}", pattern: "vcf/filtered.vcf.{gz,gz.tbi}", mode: 'move'

    input:
    path "raw.vcf"
    // setting a string value like this, ensures the input file will always have the same filename

    output:
    path "vcf/filtered.vcf.{gz,gz.tbi}"

    when:
    params.minQual > 0

    script:
    '''
    mkdir vcf
    bcftools view -i 'QUAL>10' -Oz raw.vcf > vcf/filtered.vcf.gz
    tabix vcf/filtered.vcf.gz
    '''
}

