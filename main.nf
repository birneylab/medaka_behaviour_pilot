#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Medaka behaviour pilot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

include { PREPROCESSING } from './subworkflows/preprocessing'
include { TRACKING      } from './subworkflows/tracking'
include { HMM           } from './subworkflows/hmm'

workflow  TRACK_VIDEOS {
    Channel.fromPath ( params.input_tracking )
    .splitCsv ( header: true )
    .map { [it, it.video] }
    .set { in_vid_ch }
    PREPROCESSING ( in_vid_ch )
    TRACKING ( PREPROCESSING.out )
}

workflow RUN_HMM {
    // this is decoupled from the above
    Channel.fromPath ( params.input_traj )
    .splitCsv ( header: true )
    .map { [it, it.video] }
    .first()
    .set { in_vid_ch }
    //HMM ()
}

workflow {
    // workflows are decoupled, just comment out what you don't want to run
    // and provide input files appropriately
    TRACK_VIDEOS()
    //HMM()
}