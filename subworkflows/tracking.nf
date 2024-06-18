#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Medaka behaviour pilot - tracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

process track_video {
    label "idtrackerai"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(video_in)
        )

    output:
        tuple(
            val(meta),
            path("session_${meta.id}", type: "dir"),
            // output needed so that it fails if no output, error handling not perfect in idtrackerai
            path("session_${meta.id}/trajectories/with_gaps_csv/trajectories.csv")
        )


    script:
        """
        export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib"

        { 
            idtrackerai \\
                --track \\
                --output_dir . \\
                --name ${meta.id} \\
                --video_paths ${video_in} \\
                --use_bkg ${meta.bgsub} \\
                --background_subtraction_stat ${meta.bgsub_mode} \\
                --tracking_intervals "0,${meta.video_length}" \\
                --number_of_animals 2 \\
                --intensity_ths ${meta.intensity_floor} ${meta.intensity_ceiling} \\
                --area_ths ${meta.area_floor} ${meta.area_ceiling} \\
                --protocol3_action ${meta.protocol3_action} \\
                --add_time_column_to_csv TRUE
        } || {
            # needed because idtrackerai fails when it shouldn't
            exit 0
        }
        """
}

process unpack_tracking_results {
    label "python_opencv_numpy_pandas"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(session_folder)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_traj.csv.gz"),
            path("${meta.id}_meta.json")
        )


    script:
        """
        #!/usr/bin/env python3

        import json
        import numpy as np
        import pandas as pd
        from pathlib import Path

        # sometimes the without_gaps.npy file is not created if the trajectories
        # are already without gaps
        f = Path("${session_folder}/trajectories/without_gaps.npy")
        skip_prob = False
        if not f.is_file():
            f = Path("${session_folder}/trajectories/with_gaps.npy")
            skip_prob = True # no id_probabilities if there are no gaps!
        assert f.is_file()

        # arr.item() needed to unpack 0-dimensional array
        arr = np.load(f, allow_pickle = True).item()
        traj = arr["trajectories"] # (N_frames, N_animals, 2)
        if not skip_prob:
            id_prob = arr["id_probabilities"] # (N_frames, N_animals, 1)
        else:
            id_prob = np.ones([traj.shape[0], 2, 1]) # set all id_prob to 1 if there are no gaps

        df = pd.DataFrame(
            np.concatenate(
                (
                    traj.reshape(traj.shape[0], 4),
                    id_prob.squeeze()
                ),
                axis = 1
            ),
            columns = ["x1", "y1", "x2", "y2", "id_prob1", "id_prob2"]
        )

        meta = {
            "id": "${meta.id}",
            "fps": arr["frames_per_second"],
            "version": arr["version"]
        } | arr["stats"] # | is the union operator for dicttionaries

        df.to_csv("${meta.id}_traj.csv.gz", index = False)
        with open("${meta.id}_meta.json", "w") as f_out:
            json.dump(meta, f_out)
        """
}

process assign_ref_test {
    // determine which blob is the reference fish
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(traj),
            val(cab_coords)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_traj_with_identities.csv.gz")
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")

        ref <- "${meta.ref}"
        test <- "${meta.test}"
        traj <- fread("${traj}")
        cab_coords <- "${cab_coords}"

        first_frame <- traj[!is.na(x1) & !is.na(x2) & !is.na(y1) & !is.na(y2)][1,]
        if (ref == test) {
            message("Ref is the same as test")
            stopifnot(cab_coords == "")
            idx <- 1
        } else {
            stopifnot(cab_coords %in% c("Top", "Bottom", "Left", "Right"))
            if (cab_coords == "Top") {
                idx <- which.min(c(first_frame[["y1"]], first_frame[["y2"]]))
            } else if (cab_coords == "Bottom") {
                idx <- which.max(c(first_frame[["y1"]], first_frame[["y2"]]))
            } else if (cab_coords == "Left") {
                idx <- which.min(c(first_frame[["x1"]], first_frame[["x2"]]))
            } else if (cab_coords == "Right") {
                idx <- which.max(c(first_frame[["x1"]], first_frame[["x2"]]))
            } else {
                stop("Unrecognised ref position")
            }
        }

        if (idx == 1) {
            colnames(traj)[colnames(traj) == "x1"] <- "ref_x"
            colnames(traj)[colnames(traj) == "x2"] <- "test_x"
            colnames(traj)[colnames(traj) == "y1"] <- "ref_y"
            colnames(traj)[colnames(traj) == "y2"] <- "test_y"
            colnames(traj)[colnames(traj) == "id_prob1"] <- "ref_id_prob"
            colnames(traj)[colnames(traj) == "id_prob2"] <- "test_id_prob"
        } else if (idx == 2) {
            colnames(traj)[colnames(traj) == "x1"] <- "test_x"
            colnames(traj)[colnames(traj) == "x2"] <- "ref_x"
            colnames(traj)[colnames(traj) == "y1"] <- "test_y"
            colnames(traj)[colnames(traj) == "y2"] <- "ref_y"
            colnames(traj)[colnames(traj) == "id_prob1"] <- "test_id_prob"
            colnames(traj)[colnames(traj) == "id_prob2"] <- "ref_id_prob"
        } else {
            stop("Index not recognised")
        }

        fwrite(traj, "${meta.id}_traj_with_identities.csv.gz")
        """
}

process visualise_identities {
    // check that ref and test are assigned correctly
    label "python_opencv_numpy_pandas"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(video_in),
            path(traj)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_identities.avi")
        )

    script:
        """
        #!/usr/bin/env python3

        import itertools
        import cv2 as cv
        import numpy as np
        import pandas as pd

        circle_size = 2
        font_size = 0.3
        font_width = 1
        colors = dict(
            black = (0, 0, 0), # BGR black
            red = (0, 0, 255), # BGR red
            green = (0, 255, 0), # BGR green
        )
        font = cv.FONT_HERSHEY_SIMPLEX
        of_start = ${meta.of_start}
        of_end = ${meta.of_end}
        no_start = ${meta.no_start}
        no_end = ${meta.no_end}
        adj_top_of = ${meta.adj_top_of}
        adj_right_of = ${meta.adj_right_of}
        adj_top_no = ${meta.adj_top_no}
        adj_right_no = ${meta.adj_right_no}

        cap = cv.VideoCapture("${video_in}")
        n_frames = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
        fps = int(cap.get(cv.CAP_PROP_FPS))
        w = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
        h = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
        mid_x_of = round(((w - 1) / 2) + adj_right_of)
        mid_y_of = round(((h - 1) / 2) + adj_top_of)
        mid_x_no = round(((w - 1) / 2) + adj_right_no)
        mid_y_no = round(((h - 1) / 2) + adj_top_no)

        fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')
        out = cv.VideoWriter(
            "${meta.id}_identities.avi", fourcc, fps, (w, h), isColor = True
        ) 
        
        df_dict = {"of": {}, "no": {}}
        skip_q3 = ("${video_in.baseName - "_flipped"}" == "20190616_1227_icab_kaga_R")
        for a, q in itertools.product(("of", "no"), range(1, 5)):
            f = "${video_in.baseName - "_flipped"}_{a}_q{q}_traj_with_identities.csv.gz".format(a = a, q = q)
            if skip_q3 and q == 3:
                continue
            df_dict[a][str(q)] = pd.read_csv(f)

            if a == "of":
                start = of_start
                mid_x = mid_x_of
                mid_y = mid_y_of
            elif a == "no":
                start = no_start
                mid_x = mid_x_no
                mid_y = mid_y_no
            else:
                raise AssertionError

            if q in [1, 2]:
                y_offset = 0
            elif q in [3, 4]:
                y_offset = mid_y
            else:
                raise AssertionError
            
            if q in [2, 3]:
                x_offset = 0
            elif q in [1, 4]:
                x_offset = mid_x
            else:
                raise AssertionError

            df_dict[a][str(q)].index = range(start, df_dict[a][str(q)].shape[0] + start)
            df_dict[a][str(q)][["ref_x", "test_x"]] += x_offset
            df_dict[a][str(q)][["ref_y", "test_y"]] += y_offset

        def add_coloured_split_lines(frame, i):
            # Add colored split lines that indicate the part of assay we are in
            if i < of_start:
                line_color = colors["red"]
                mid_x = mid_x_of
                mid_y = mid_y_of
            elif i < of_end:
                line_color = colors["green"]
                mid_x = mid_x_of
                mid_y = mid_y_of
            elif i < no_start:
                line_color = colors["red"]
                mid_x = mid_x_no
                mid_y = mid_y_no
            elif i < no_end:
                line_color = colors["green"]
                mid_x = mid_x_no
                mid_y = mid_y_no
            else:
                line_color = colors["red"]
                mid_x = mid_x_no
                mid_y = mid_y_no
            start_point = (mid_x, 0)
            end_point = (mid_x, h)
            thickness = 1
            frame = cv.line(frame, start_point, end_point, line_color, thickness)
            start_point = (0, mid_y)
            end_point = (w, mid_y)
            frame = cv.line(frame, start_point, end_point, line_color, thickness)

            return frame

        def add_fish_coords(frame, i):
            if i < of_start:
                return frame
            elif i < of_end:
                a = "of"
            elif i < no_start:
                return frame
            elif i < no_end:
                a = "no"
            else:
                return frame

            for q in range(1, 5):
                if skip_q3 and q == 3:
                    continue
                if not df_dict[a][str(q)].loc[i][["ref_x", "ref_y"]].isna().any():
                    ref_x = int(df_dict[a][str(q)].loc[i]["ref_x"])
                    ref_y = int(df_dict[a][str(q)].loc[i]["ref_y"])
                    frame = cv.circle(
                        frame, (ref_x, ref_y), circle_size, colors["red"], -1
                    )
                    frame = cv.putText(
                        frame,
                        "ref",
                        (ref_x, ref_y),
                        font,
                        font_size,
                        colors["red"],
                        font_width,
                    )

                if not df_dict[a][str(q)].loc[i][["test_x", "test_y"]].isna().any():
                    test_x = int(df_dict[a][str(q)].loc[i]["test_x"])
                    test_y = int(df_dict[a][str(q)].loc[i]["test_y"])
                    frame = cv.circle(
                        frame, (test_x, test_y), circle_size, colors["black"], -1
                    )
                    
                    frame = cv.putText(
                        frame,
                        "test",
                        (test_x, test_y),
                        font,
                        font_size,
                        colors["black"],
                        font_width,
                    )

            return frame

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
        
        for i in range(n_frames):
            ret, frame = cap.read()
            assert ret
            frame = add_coloured_split_lines(frame, i)
            frame = add_fish_coords(frame, i)
            frame = add_frame_counter(frame, i)
            out.write(frame)
            
        cap.release()
        out.release()
        """
}

process aggregate_tracking_stats {
    // creates a dataset of tracking statistics for further exploration
    // estimated_accuracy: ?
    // estimated_accuracy_after_interpolation: mean id probability ignoring NA
    // percentage_identified: fraction of frames*animals with non-NA coordinates
    // estimated_accuracy_identified: mean id probability ignoring NA for frames*animals with non-NA coordinates
    label "python_opencv_numpy_pandas"

    input:
        path(stats)

    output:
        path("tracking_stats.csv.gz")

    script:
        """
        #!/usr/bin/env python3

        import json
        import glob
        import pandas as pd

        f_list = glob.glob("*.json")
        df = pd.DataFrame([json.load(open(f)) for f in f_list])
        df.to_csv("tracking_stats.csv.gz", index = False)
        """
}


workflow TRACKING {
    take:
        split_vids
        raw_vids
    
    main:
        Channel.fromPath ( params.idtrackerai_params )
        .splitCsv ( header: true )
        .map { meta -> [meta.id, meta] }
        .set { tracking_params }
        
        split_vids
        .map { meta, vid -> [meta.id, meta, vid] }
        .join ( tracking_params, by: 0, failOnDuplicate: true, failOnMismatch: true )
        .map { key, meta, vid, meta2 -> [meta + meta2, vid] }
        .set { track_video_in_ch }

        track_video ( track_video_in_ch )
        
        track_video.out
        .map { meta, session, traj -> [meta, session] }
        .set { tracking_sessions }

        unpack_tracking_results ( tracking_sessions )
        
        Channel.fromPath ( params.cab_coords )
        .splitCsv ( header: true )
        .map { [it.id, it.cab_coords] }
        .set { cab_coords }
        
        unpack_tracking_results.out
        .map { meta, traj, stats -> [meta.id, meta, traj] }
        .join ( cab_coords, by: 0, failOnDuplicate: true, failOnMismatch: true )
        .map { key, meta, traj, cab_coords -> [meta, traj, cab_coords] }
        .set { assign_ref_test_in }

        assign_ref_test ( assign_ref_test_in )
        aggregate_tracking_stats ( unpack_tracking_results.out.map { it[2] }.collect() )

        assign_ref_test.out
        .map {
             meta, traj ->
             [ meta.id.replaceFirst(/_(of|no)_q(1|2|3|4)$/, ""), traj ]
        }
        .groupTuple ( by: 0 )
        .join ( raw_vids.map { meta, vid -> [ meta.id, meta, vid ] } , by: 0 )
        .map { key, traj, meta, vid -> [ meta, vid, traj ] }
        .set { visualise_identities_in }
        visualise_identities ( visualise_identities_in )
}