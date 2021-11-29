#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         =============================================
          E C S E Q - S N P c a l l   P I P E L I N E
         =============================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run ecseq/dnaseq [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Provide the directory containing BAM file(s) in "*.bam" format

              --reference [path/to/ref.fa]    [REQUIRED] Provide the path to the reference genome in fasta format

              --output [STR]                  A string that can be given to name the output directory. [default: "."]


         Options: MODIFIERS
              --markDups                      Mark PCR duplicates with Picard MarkDuplicates. [default: off]

              --bamQC                         Generate bamQC report of alignments. [default: off]


         Options: FILTERING
              --minQual                       Minimum variant quality threshold. [default: 0]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run ecseq/SNPcall \
              --input /path/to/input/dir \
              --reference /path/to/genome.fa \
              --markDups --bamQC --minQual 20

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         =============================================
          E C S E Q - S N P c a l l   P I P E L I N E
         =============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
fasta = file("${params.reference}", checkIfExists: true, glob: false)
fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
bam_path = "${params.input}/*.bam"



// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ============================================="
log.info "          E C S E Q - S N P c a l l   P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         =============================================" }
else {
log.info "         =============================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${workflow.profile.tokenize(",").contains("test") ? "-" : "${bam_path}"}"
log.info "         reference    : ${params.reference}"
log.info "         output dir   : ${params.output}"
log.info "         mode         : ${params.minQual > 0 ? "Filtered" : "RAW"}"
log.info "         QC options   : ${params.markDups ? "markDups " : ""}${params.bamQC ? "bamQC" : ""}"
log.info ""
log.info "         ============================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, i.e. the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './lib/functions.nf' params(bamPaths: params.bamPaths)
        BAMS = check_test_data(params.bamPaths)

} else {

    // STAGE READS CHANNELS # this defines the normal input when test profile is not in use
    BAMS = Channel
        .fromPath(bam_path)
        .ifEmpty{ exit 1, "ERROR: cannot find valid read files in dir: ${params.input}\n \
        The pipeline will expect BAM files in *.bam format"}
        .take(params.take.toInteger())

}



////////////////////
// BEGIN PIPELINE //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'DNAseq'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */

// INCLUDES # here you must give the relevant process files from the lib directory 
include {samtools_merge;Picard_MarkDuplicates;Freebayes;bcftools;bamQC} from './lib/process.nf' params(params)

// SUB-WORKFLOWS
workflow 'SNPcall' {

    // take the initial Channels and paths
    take:
        BAMS
        fasta
        fai

    // here we define the structure of our workflow i.e. how the different processes lead into each other
    // eg. process(input1, input2, etc.)
    // eg. process.out[0], process.out[1], etc.
    // index numbers [0],[1],etc. refer to different outputs defined for processes in process.nf
    // ALWAYS PAY ATTENTION TO CARDINALITY!!
    main:

	// we should use the collect operator to get all BAM files into a single item emitted by the Channel
        samtools_merge(BAMS.collect())

	// Picard MarkDuplicates only ever takes the output Channel from samtools_merge
	Picard_MarkDuplicates(samtools_merge.out)

	// Freebayes and bamQC run differently depending on whether or not we use --markDups
	if(params.markDups){
	Freebayes(Picard_MarkDuplicates.out[0],fasta,fai)
	bamQC(Picard_MarkDuplicates.out[0])
	} else {
	Freebayes(samtools_merge.out,fasta,fai)
	bamQC(samtools_merge.out)
	}

	// bcftools filtering only ever runs on the output from Freebayes
	bcftools(Freebayes.out)
}

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        SNPcall(BAMS, fasta, fai)

}


//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
