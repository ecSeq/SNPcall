# ecSeq-SNPcall Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Preprocessing](#preprocessing) - merging individual bam files into one
* [Quality Control](#quality-control) - marking PCR duplicates and/or generating bamQC reports
* [Variant calling](#variant-calling) - joint variant calling with Freebayes
* [Pipeline Info](#pipeline-info) - reports from nextflow about the pipeline run

## Preprocessing
Input bam files will be merged using samtools merge, which differentiates each set of alignments by Read Group.

## Quality Control
Following preprocessing, the pipeline will mark PCR duplicates in the merged bam file. The pipeline will generate a QC report for either the merged or marked bam file.

**Output directory: `./bam`**

## Variant calling
Depending on which options are specified to the pipeline, either the marked or merged bam file will be taken for joint variant calling with Freebayes to produce a `raw.vcf` file. Optional quality filtering is performed with bcftools to produce a `filtered.vcf.gz` file.

**Output directory: `./vcf`**

## Pipeline Info
Nextflow has several built-in reporting tools that give information about the pipeline run.

**Output directory: `./`**

* `dag.svg`
  * DAG graph giving a diagrammatic view of the pipeline run.
  * NB: If [Graphviz](http://www.graphviz.org/) was not installed when running the pipeline, this file will be in [DOT format](http://www.graphviz.org/content/dot-language) instead of SVG.
* `report.html`
  * Nextflow report describing parameters, computational resource usage and task bash commands used.
* `timeline.html`
  * A waterfall timeline plot showing the running times of the workflow tasks.
* `trace.txt`
  * A text file with machine-readable statistics about every task executed in the pipeline.
