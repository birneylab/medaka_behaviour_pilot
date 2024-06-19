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
    // compute metrics to feed to the hmm from the trajectories:
    // angles and distances
    label "r_tidyverse_datatable"

    input:
        tuple(
            val(meta),
            path(traj),
            val(time_step)
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

        get_dist <- function(x1, x2, y1, y2) {
            sqrt((x1 - x2)^2 + (y1 - y2)^2)
        }

        get_angle <- function(x, xlag1, xlag2, y, ylag1, ylag2) {
            # segments defining the angle
            x_bc <- x - xlag1
            y_bc <- y - ylag1
            x_ab <- xlag1 - xlag2
            y_ab <- ylag1 - ylag2
            # dot product is the element-wise product
            dot <- (x_ab * x_bc) + (y_ab * y_bc)
            # determinant is the difference of the diagonals
            det <- (x_ab * y_bc) - (y_ab * x_bc)
            # det is proportional to sin, dot is proportional to cos, with the same constant
            rad <- atan2(det, dot)
            return(rad)
        }

        step_size_s <- ${time_step}
        step_size_frames <- round(step_size_s * ${meta.fps})
        
        df <- fread("${traj}")
        df[, frame_n := .N]
        df[, time_s := frame_n/${meta.fps}]
        df[ref_x_lag := lag(ref_x, step_size_frames)]
        df[ref_y_lag := lag(ref_y, step_size_frames)]
        df[test_x_lag := lag(test_x, step_size_frames)]
        df[test_y_lag := lag(test_y, step_size_frames)]
        df[ref_x_lag2 := lag(ref_x, step_size_frames * 2)]
        df[ref_y_lag2 := lag(ref_y, step_size_frames * 2)]
        df[test_x_lag2 := lag(test_x, step_size_frames * 2)]
        df[test_y_lag2 := lag(test_y, step_size_frames * 2)]
        df[ref_distance := get_dist(ref_x, ref_x_lag, ref_y, ref_y_lag)]
        df[test_distance := get_dist(test_x, test_x_lag, test_y, test_y_lag)]
        df[ref_angle := get_angle(ref_x, ref_x_lag, ref_x_lag2, ref_y, ref_y_lag, ref_y_lag2)]
        df[test_angle := get_angle(test_x, test_x_lag, test_x_lag2, test_y, test_y_lag, test_y_lag2)]
        out <- df[
            , .(
                frame_n, time_s, ref_x, ref_y, test_x, test_y, ref_distance, test_distance, ref_angle, test_angle
            )
        ]

        fwrite(df, "${meta.id}_metrics.csv.gz")
        """
}

workflow HMM {
    take:
        in_ch
    
    main:
        in_ch.view()
        //compute_metrics ( in_ch )
}