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
            path("${meta.id}_hmm.csv.gz")
        )

    script:
        """
        #!/usr/bin/env python3

        from hmmlearn import hmm
        import pandas as pd
        import numpy as np
        import glob
        import re
        import random

        random.seed(1)
        hmm_seed = random.randint(0, int(1e6))

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
            df_ref = df[["id", "frame_n", "time_s", "ref_distance", "ref_angle"]].rename(columns = col_renamer).dropna()
            df_test = df[["id", "frame_n", "time_s", "test_distance", "test_angle"]].rename(columns = col_renamer).dropna()
            df_ref["id"] += "_ref"
            df_test["id"] += "_test"

            return pd.concat([df_ref, df_test], ignore_index = True)

        f_list = glob.glob("*_metrics.csv.gz")
        df = pd.concat(
            map(read_metric, f_list),
            ignore_index = True
        ).sort_values(
            by = ["id", "frame_n"]
        )

        X = df[["distance", "angle"]].to_numpy()
        l = df.groupby("id").size().to_numpy()
        assert l.sum() == X.shape[0]

        model = hmm.GaussianHMM(n_components=${n_states}, covariance_type="diag", n_iter=100, random_state = hmm_seed)
        model.fit(X, lengths = l)

        out = df
        out["hmm_state"] = model.predict(X, lengths = l)
        out.to_csv("${meta.id}_hmm.csv.gz", index = False)
        """
}

process hmm_cross_validation {
    label "python_hmmlearn_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(metrics),
            val(n_states),
            path(cv_splits)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_hmm_cross_validation.csv.gz")
        )

    script:
        """
        #!/usr/bin/env python3

        from hmmlearn import hmm
        import pandas as pd
        import numpy as np
        import glob
        import re
        import random

        random.seed(1)

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
            df_ref = df[["id", "frame_n", "time_s", "ref_distance", "ref_angle"]].rename(columns = col_renamer).dropna()
            df_test = df[["id", "frame_n", "time_s", "test_distance", "test_angle"]].rename(columns = col_renamer).dropna()
            df_ref["id"] += "_ref"
            df_test["id"] += "_test"

            return pd.concat([df_ref, df_test], ignore_index = True)

        def train_hmm(f_df, cv_class):
            hmm_seed = random.randint(0, int(1e6))
            
            f_list = f_df.loc[f_df["cv_fold"] == cv_class, "f"]
            df = pd.concat(
                map(read_metric, f_list),
                ignore_index = True
            ).sort_values(
                by = ["id", "frame_n"]
            )
            df["cv_fold"] = cv_class

            X = df[["distance", "angle"]].to_numpy()
            l = df.groupby("id").size().to_numpy()
            assert l.sum() == X.shape[0]

            model = hmm.GaussianHMM(
                n_components=${n_states}, covariance_type="diag", n_iter=1, random_state = hmm_seed
            )
            model.fit(X, lengths = l)

            ret = {
                "model": model,
                "X": X,
                "l": l,
                "df": df
            }

            return ret

        def hmm_predict(d, self, val):
            df = d[self]["df"]
            X = d[self]["X"]
            l = d[self]["l"]

            df["hmm_state_self"] = list(map(
                lambda s: self + "_" + str(s),
                d[self]["model"].predict(X, lengths = l)
            ))
            df["hmm_state_val"] = list(map(
                lambda s: val + "_" + str(s),
                d[val]["model"].predict(X, lengths = l)
            ))

            return df

        cv_splits = pd.read_csv("${cv_splits}").set_index("id")
        f_list_full = glob.glob("*_metrics.csv.gz")
        f_df = pd.DataFrame(
            {
                "id": map(lambda f: re.sub("_metrics.csv.gz", "", f), f_list_full),
                "f": f_list_full,
            }
        ).join(
            cv_splits,
            on = "id",
            validate = "one_to_one"
        )
        
        res = {
            "A": train_hmm(f_df, "A"),
            "B": train_hmm(f_df, "B")
        }
        out = pd.concat(
            [
                hmm_predict(res, self = "A", val = "B"),
                hmm_predict(res, self = "B", val = "A")
            ]
        )
        out.to_csv("${meta.id}_hmm_cross_validation.csv.gz", index = False)
        """
}

process hmm_concordance {
    label "r_tidyverse_datatable"

    input:
        tuple(
            val(meta),
            path(hmm_cv)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_concordance.csv.gz")
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")
        
        df <- fread("${hmm_cv}")
        concordance <- list()
        
        # process separately the 2 folds
        for (the_cv_fold in c("A", "B")) {
            tmp <- df[cv_fold == the_cv_fold]
            # order self states from most to least populated
            self_state_counts <- tmp[, .(n = .N), by = hmm_state_self][order(n, decreasing = TRUE)]
            state_match <- data.table(
                self_state = self_state_counts[["hmm_state_self"]],
                val_state = "NA", # must be a character!
                cv_fold = the_cv_fold
            )
            
            # get the state matches
            for (the_self_state in state_match[["self_state"]]) {
                # get the val_state with the highest overlap to self_state
                val_state_counts <- tmp[
                    hmm_state_self == the_self_state & !(hmm_state_val %in% state_match[["val_state"]]),
                    .(n = .N),
                    by = hmm_state_val
                ]
                
                if (nrow(val_state_counts) > 0) {
                    the_val_state <- val_state_counts[
                        order(n, decreasing = TRUE)[1], hmm_state_val
                    ]
                } else {
                    the_val_state <- tmp[
                        !(hmm_state_val %in% state_match[["val_state"]]),
                        unique(hmm_state_val)[1]
                    ]
                }

                # assign the match
                state_match[self_state == the_self_state, val_state := ..the_val_state]
            }

            stopifnot(sort(state_match[["val_state"]]) == sort(unique(tmp[["hmm_state_val"]])))
            stopifnot(sort(state_match[["self_state"]]) == sort(unique(tmp[["hmm_state_self"]])))

            # recode the states
            tmp[, hmm_state_self := match(hmm_state_self, ..state_match[["self_state"]])]
            tmp[, hmm_state_val := match(hmm_state_val, ..state_match[["val_state"]])]

            # compute concordance
            concordance[[the_cv_fold]] <- tmp[, mean(hmm_state_self == hmm_state_val)]
        }

        out <- data.table(
            time_step = ${meta.time_step},
            n_states = ${meta.n_states},
            concordance_A = concordance[["A"]],
            concordance_B = concordance[["B"]]
        )
        fwrite(out, "${meta.id}_concordance.csv.gz")
        """
}

process combine_concordance {
    label "r_tidyverse_datatable"

    input:
        path(concordance)

    output:
        path("hmm_concordance.csv.gz")

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")
        
        f_list <- list.files(pattern = "*_concordance.csv.gz")
        df <- lapply(f_list, fread) |> rbindlist()
        fwrite(df, "concordance.csv.gz")
        """
}


workflow HMM {
    take:
        traj
        split_vids
    
    main:
        traj.combine ( params.time_step )
        .filter { it[2] == 0.08 }
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
        .filter { it[2] == 14 }
        .map {
            meta, metrics, n_states ->
            def new_meta = [
                id: "time_step${meta.time_step}_n_states${n_states}",
                time_step: meta.time_step,
                n_states: n_states
            ]
            [ new_meta, metrics, n_states ]
        }
        .groupTuple ( by: [0, 2] )
        .set { hmm_in }
        hmm_cross_validation ( hmm_in.combine ( [ params.hmm_cv_splits ] ) )
        hmm_concordance ( hmm_cross_validation.out )
        combine_concordance ( hmm_concordance.out.map { meta, f -> f }.collect() )

        hmm_in.filter {
            meta, metrics, n_states ->
            meta.n_states == params.n_states_selected && meta.time_step == params.time_step_selected
        }
        .set { hmm_in_selected }
        //run_hmm ( hmm_in_selected )
}