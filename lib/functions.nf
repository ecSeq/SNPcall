#!/usr/bin/env nextflow
// This file for loading custom functions into the main.nf script (separated for portability)

// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(bamPaths) {

    // Set BAMS testdata
    BAMS = Channel.from(bamPaths)
                  .ifEmpty { exit 1, "test profile bamPaths was empty - no input files supplied" }

    // Return BAMS channel
    return BAMS
}
