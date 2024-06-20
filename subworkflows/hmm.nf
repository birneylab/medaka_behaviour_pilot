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

        get_angle <- function(xlag, x, xlead, ylag, y, ylead) {
            # this is equivalent to the difference in heading between the incoming and outgoing segments
            # segments defining the angle
            x_bc <- xlead - x
            y_bc <- ylead - y
            x_ab <- x - xlag
            y_ab <- y - ylag
            # dot product is the element-wise product
            dot <- (x_ab * x_bc) + (y_ab * y_bc)
            # determinant is the difference of the diagonals
            det <- (x_ab * y_bc) - (y_ab * x_bc)
            # det is proportional to sin, dot is proportional to cos, with the same constant
            rad <- atan2(det, dot)
            return(rad)
        }

        get_heading <- function(xlag, x, ylag, y) {
            px <- x - xlag
            py <- y - ylag
            rad <- atan2(py, px)
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
        df[, ref_x_lead := shift(ref_x, -step_size_frames)]
        df[, ref_y_lead := shift(ref_y, -step_size_frames)]
        df[, test_x_lead := shift(test_x, -step_size_frames)]
        df[, test_y_lead := shift(test_y, -step_size_frames)]
        df[, ref_distance := get_dist(ref_x, ref_x_lag, ref_y, ref_y_lag)]
        df[, test_distance := get_dist(test_x, test_x_lag, test_y, test_y_lag)]
        df[, ref_angle := get_angle(ref_x_lag, ref_x, ref_x_lead, ref_y_lag, ref_y, ref_y_lead)]
        df[, test_angle := get_angle(test_x_lag, test_x, test_x_lead, test_y_lag, test_y, test_y_lead)]
        df[, ref_heading := get_heading(ref_x_lag, ref_x, ref_y_lag, ref_y)]
        df[, test_heading := get_heading(test_x_lag, test_x, test_y_lag, test_y)]
        out <- df[
            , .(
                frame_n, time_s, ref_x, ref_y, test_x, test_y, ref_distance, test_distance, ref_angle, test_angle, ref_heading, test_heading
            )
        ]

        fwrite(df, "${meta.id}_metrics.csv.gz")
        """
}

process visualise_metrics {
    label "python_opencv_numpy_pandas"
    queue "datamover"

    input:
        tuple(
            val(meta),
            path(metrics),
            path(video_in)
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
        fps = int(cap.get(cv.CAP_PROP_FPS))
        w = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
        h = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
        
        fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')
        out = cv.VideoWriter(
            "${meta.id}_metrics.avi", fourcc, fps, (w, h), isColor = True
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
                x1 = df.iloc[i][["{}_x".format(fish), "{}_y".format(fish)]].to_numpy()
                theta = df.iloc[i]["{}_angle".format(fish)]
                heading = df.iloc[i]["{}_heading".format(fish)]
                delta = df.iloc[i]["{}_distance".format(fish)]
                # I want to show the turn with respect to the previous heading, not from the x axis
                dx_dy = np.array([np.cos(theta + heading), np.sin(theta + heading)]) * delta
                x2 = x1 + (dx_dy * scale_factor)

                if np.isnan(x1).any() or np.isnan(x2).any():
                    continue

                pt1 = x1.round().astype(int)
                pt2 = x2.round().astype(int)

                frame = cv.arrowedLine(
                    frame,
                    pt1,
                    pt2,
                    colors[color_name],
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

process run_hmm {
    label "python_hmmlearn_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(metrics),
            val(n_states)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_states${meta.n_states}_step${meta.time_step}_hmm.csv.gz")
        )

    script:
        """
        #!/usr/bin/env python3

        from hmmlearn import hmm
        import pandas as pd
        import numpy as np
        import glob
        import re

        col_renamer = {
            "ref_distance": "distance",
            "test_distance": "distance",
            "ref_angle": "angle",
            "test_angle": "angle",
        }

        def read_metric(f):
            id = re.sub("_metrics.csv.gz", "", f)
            df = pd.read_csv(f)
            df["id"] = id
            df_ref = df[["id", "ref_distance", "ref_angle"]].rename(columns = col_renamer)
            df_test = df[["id", "test_distance", "test_angle"]].rename(columns = col_renamer)
            df_ref["id"] += "ref"
            df_test["id"] += "test"

            return pd.concat([df_ref, df_test])

        f_list = glob.glob("*_metrics.csv.gz")
        df = pd.concat(
            map(read_metric, f_list)
        )
        X = df[["distance", "angle"]].to_numpy()
        l = df.groupby("id").size().to_numpy()

        model = hmm.GaussianHMM(n_components=${n_states}, covariance_type="diag", n_iter=100)
        model.fit(X, lengths = l)

        out = df[["id"]]
        out[["hmm_state"]] = model.predict(X, lengths = l)
        out.to_csv("${meta.id}_states${meta.n_states}_step${meta.time_step}_hmm.csv.gz", index = False)
        """
}


workflow HMM {
    take:
        traj
        split_vids
    
    main:
        traj.combine ( params.time_interval )
        .map{
            meta, traj, time_step ->
            def new_meta = meta.clone()
            new_meta.time_step = time_step
            [ new_meta, traj, time_step ]
        }
        .set { metrics_in }
        compute_metrics ( metrics_in )

        compute_metrics.out
        .map { meta, metrics -> [ meta.id, meta, metrics ] }
        .combine( split_vids, by: 0 )
        .map { id, meta, metrics, vid -> [ meta, metrics, vid ] }
        .set { visualise_metrics_in }
        visualise_metrics ( visualise_metrics_in )
        
        compute_metrics.out
        .combine( params.n_states )
        .map {
            meta, metrics, n_states ->
            def new_meta = meta.clone()
            new_meta.n_states = n_states
            def group_key = [meta.time_step, n_states]
            [ group_key, new_meta, metrics, n_states ]
        }
        .groupTuple(by: 0)
        .map { group_key, meta, metrics, n_states -> [ meta, metrics, n_states ] }
        .set { hmm_in }
        //run_hmm ( hmm_in )
}