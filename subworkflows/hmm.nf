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
    queue "datamover"

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
        df[, frame_n := 1:.N]
        df[, time_s := frame_n/${meta.fps}]
        df[, ref_x_lag := shift(ref_x, step_size_frames)]
        df[, ref_y_lag := shift(ref_y, step_size_frames)]
        df[, test_x_lag := shift(test_x, step_size_frames)]
        df[, test_y_lag := shift(test_y, step_size_frames)]
        df[, ref_x_lag2 := shift(ref_x, step_size_frames * 2)]
        df[, ref_y_lag2 := shift(ref_y, step_size_frames * 2)]
        df[, test_x_lag2 := shift(test_x, step_size_frames * 2)]
        df[, test_y_lag2 := shift(test_y, step_size_frames * 2)]
        df[, ref_distance := get_dist(ref_x, ref_x_lag, ref_y, ref_y_lag)]
        df[, test_distance := get_dist(test_x, test_x_lag, test_y, test_y_lag)]
        df[, ref_angle := get_angle(ref_x, ref_x_lag, ref_x_lag2, ref_y, ref_y_lag, ref_y_lag2)]
        df[, test_angle := get_angle(test_x, test_x_lag, test_x_lag2, test_y, test_y_lag, test_y_lag2)]
        out <- df[
            , .(
                frame_n, time_s, ref_x, ref_y, test_x, test_y, ref_distance, test_distance, ref_angle, test_angle
            )
        ]

        fwrite(df, "${meta.id}_metrics.csv.gz")
        """
}

process visualise_metrics {
    label "python_opencv_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(metrics)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_metrics.avi")
        )

    script:
        """
        #!/usr/bin/env python3

        import cv2 as cv
        import numpy as np
        import pandas as pd
        
        cap = cv.VideoCapture("${video_in}")
        n_frames = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
        w = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
        h = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
        
        fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')
        out = cv.VideoWriter(
            "${meta.id}_identities.avi", fourcc, fps, (w, h), isColor = True
        ) 

        colors = dict(
            black = (0, 0, 0), # BGR black
            red = (0, 0, 255), # BGR red
        )
        font = cv.FONT_HERSHEY_SIMPLEX
        font_size = 0.3
        font_width = 1
        scale_factor = 10

        df = pd.read_csv("${metrics}")

        def add_frame_counter(frame, i):
            frame = cv.putText(
                frame,
                str(i),
                (0, 10),
                font,
                font_size,
                colors["black"],
                font_width,
            )

            return frame

        def add_metrics(frame, i):
            for fish, color_name in zip(["test", "ref"], ["black", "red"]):
                pt1 = df.iloc[i][["{n}_x".format(fish), "{n}_y".format(fish)]].to_numpy()
                theta = df.iloc[i]["{n}_angle".format(fish)]
                delta = df.iloc[i]["{n}_distance".format(fish)]
                dx_dy = np.array([np.cos(theta), np.sin(theta)]) * delta
                pt2 = pt1 + (dx_dy * scale_factor)

                frame = cv.arrowedLine(
                    frame,
                    pt1,
                    pt2,
                    colors[[color_name]],
                )

            return frame
        
        for i in range(n_frames):
            ret, frame = cap.read()
            assert ret
            frame = add_frame_counter(frame, i)
            frame = add_metrics(frame, i)
            out.write(frame)
            
        cap.release()
        out.release()
        """
}


workflow HMM {
    take:
        in_ch
    
    main:
        in_ch.combine ( params.time_interval )
        .map{
            meta, traj, time_step ->
            def new_meta = meta.clone()
            new_meta.time_step = time_step
            [ new_meta, traj, time_step ]
        }
        .filter { meta.id == "20190611_1331_icab_icab_R_no_q1" }
        .set { metrics_in }
        compute_metrics ( metrics_in )
}