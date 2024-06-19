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
    TRACKING ( PREPROCESSING.out.split_videos, PREPROCESSING.out.raw_videos )
}

workflow RUN_HMM {
    // this is voluntarily decoupled from the tracking
    Channel.fromPath ( params.input_hmm )
    .splitCsv ( header: true )
    .map { [ it.id, it ] }
    .set { traj }
    
    Channel.fromPath ( params.tracking_stats )
    .splitCsv ( header: true )
    .map { [ it.id, it ] }
    .set { stats }
    
    traj.join ( stats, by: 0, failOnDuplicate: true, failOnMismatch: true )
    .map {
        key, traj_map, stat_map ->
        def meta = [
            id: key,
            fps: stat_map.fps,
        ]
        [ meta, traj_map.traj_file ] 
    }
    .set { hmm_in }

    HMM ( hmm_in )
}

workflow {
    // workflows are decoupled, just comment out what you don't want to run
    // and provide input files appropriately
    //TRACK_VIDEOS()
    RUN_HMM()
}