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
        
        df <- fread("${traj}")
        intervals <- seq(0, nrow(df), ${meta.fps} * ${params.time_interval})
        # non-overlapping intervals instead than rolling window because this is how Ian did before
        df <- df[intervals]
        df[ref_x_lag := lag(ref_x)]
        df[ref_y_lag := lag(ref_y)]
        df[test_x_lag := lag(test_x)]
        df[test_y_lag := lag(test_y)]
        df[ref_x_lag2 := lag(ref_x, 2)]
        df[ref_y_lag2 := lag(ref_y, 2)]
        df[test_x_lag2 := lag(test_x, 2)]
        df[test_y_lag2 := lag(test_y, 2)]
        df[ref_distance := get_dist(ref_x, ref_x_lag, ref_y, ref_y_lag)]
        df[test_distance := get_dist(test_x, test_x_lag, test_y, test_y_lag)]
        df[ref_angle := get_angle(ref_x, ref_x_lag, ref_x_lag2, ref_y, ref_y_lag, ref_y_lag2)]
        df[test_angle := get_angle(test_x, test_x_lag, test_x_lag2, test_y, test_y_lag, test_y_lag2)]

        fwrite(df, "${meta.id}_metrics.csv.gz")
        """
}

workflow HMM {
    take:
        in_ch
    
    main:
        in_ch.view()
}