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
            path(traj)
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
        cab_coords <- ${meta.ref == meta.test ? "NA" : "\"${meta.cab_coords}\"" }

        first_frame <- traj[!is.na(x1) & !is.na(x2) & !is.na(y1) & !is.na(y2)][1,]
        if (ref == test) {
            message("Ref is the same as test")
            idx <- 1
        } else {
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
            path("${meta.id}_identities.png")
        )

    script:
        """
        #!/usr/bin/env python3

        import cv2 as cv
        import numpy as np
        import pandas as pd

        circle_size = 2
        font_size = 0.3
        font_width = 1
        colors = dict(
            ref = (0, 0, 255), # BGR red
            test = (0, 0, 0), # BGR black
        )
        font = cv.FONT_HERSHEY_SIMPLEX
        cap = cv.VideoCapture("${video_in}")
        df = pd.read_csv("${traj}")
        df["frame"] = range(0, df.shape[0])
        df = df.dropna()
        frame_number = int(df.iloc[0]["frame"])
        cap.set(cv.CAP_PROP_POS_FRAMES, frame_number)
        ret, frame = cap.read()
        ref_x = int(df.iloc[0]["ref_x"])
        ref_y = int(df.iloc[0]["ref_y"])
        test_x = int(df.iloc[0]["test_x"])
        test_y = int(df.iloc[0]["test_y"])

        # add points showing the test and reference fish
        cv.circle(
            frame, (ref_x, ref_y), circle_size, colors["ref"], -1
        )
        cv.circle(
            frame, (test_x, test_y), circle_size, colors["test"], -1
        )
        cv.putText(
            frame,
            "ref",
            (ref_x, ref_y),
            font,
            font_size,
            colors["ref"],
            font_width,
        )
        cv.putText(
            frame,
            "test",
            (test_x, test_y),
            font,
            font_size,
            colors["test"],
            font_width,
        )

        # Write frame
        cv.imwrite("${meta.id}_identities.png", frame)
        cap.release()
        """
}

process aggregate_tracking_stats {
    // creates a dataset of tracking statistics for further exploration
    label "python_opencv_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(stats)
        )

    output:
        path("tracking_stats.csv.gz")

    script:
        """
        #!/usr/bin/env python3

        import json
        import glob

        f_list = glob.glob("*.json")
        df = pd.DataFrame([json.load(f) for f in f_list])
        df.to_csv("tracking_stats.csv.gz", header = False)
        """
}


workflow TRACKING {
    take:
        split_vids
    
    main:
        Channel.fromPath ( params.idtrackerai_params )
        .splitCsv ( header: true )
        .map { meta -> [meta.id, meta] }
        .set { tracking_params }
        split_vids
        .map { meta, vid -> [meta.id, meta, vid] }
        .join ( tracking_params, by: 0 )
        .map { key, meta, vid, meta2 -> [meta + meta2, vid] }
        .set { track_video_in_ch }
        track_video ( track_video_in_ch )
        track_video.out
        .map { meta, session, traj -> [meta, session] }
        .set { tracking_sessions }
        unpack_tracking_results ( tracking_sessions )
        assign_ref_test ( unpack_tracking_results.out.map { it[0,1] } )
        visualise_identities ( track_video_in_ch.join ( assign_ref_test.out, by: 0 ) )
        aggregate_tracking_stats ( unpack_tracking_results.out.map { it[2] }.collect() )
}