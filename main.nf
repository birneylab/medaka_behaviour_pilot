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

workflow  {
    Channel.fromPath ( params.input )
    .splitCsv ( header: true )
    .map { [it, it.video] }
    .first()
    .set { in_ch }
    PREPROCESSING ( in_ch )
}