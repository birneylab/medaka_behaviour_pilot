#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Medaka behaviour pilot - hmm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

process compute_metrics {
    // compute metrics to feed to the hmm from the trajectories
    label "r_tidyverse_datatable"

    input:
        tuple(
            val(meta),
            path(traj)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_metrics.csv.gz")
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")
        
        df <- fread("${traj}")
        interval_frames <- ${meta.fps} * ${params.time_interval}
        """
}

workflow HMM {
    take:
        in_ch
    
    main:
        in_ch.view()
}