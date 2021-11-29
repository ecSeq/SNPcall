[<img width="200" align="right" src="docs/images/ecseq.jpg">](https://www.ecseq.com)
[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/ecseq/snpcall.svg)](https://hub.docker.com/r/ecseq/snpcall)

ecSeq-SNPcall Pipeline
======================

**ecSeq/SNPcall** is a simple bioinformatics analysis pipeline for merging DNAseq alignments and performing joint variant calling.

The workflow processes a collection of bam files using [samtools](https://github.com/samtools/samtools), producing a merged file which can optionally be taken forward to have PCR duplicates marked with [Picard MarkDuplicates](https://broadinstitute.github.io/picard/) and/or produce QC reports with [Qualimap bamQC](http://qualimap.conesalab.org/). Joint variant calling is performed with [Freebayes](https://github.com/freebayes/freebayes), with optional basic quality filtering using [bcftools](https://github.com/samtools/bcftools). 

> See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://www.nextflow.io/)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run ecseq/snpcall -profile test,<docker|singularity|conda>
```

iv. Start running your own analysis!

```bash
nextflow run ecseq/snpcall -profile <docker|singularity|conda> \
--input /path/to/bam/dir --reference /path/to/genome.fa \
<--markDups|--bamQC|--minQual 20>
```

> See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.


### Credits

These scripts were originally written for use by [ecSeq Bioinformatics GmbH](https://www.ecseq.com), by Adam Nunn ([@bio15anu](https://github.com/bio15anu)).
