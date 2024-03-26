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

    input:
        tuple(
            val(meta),
            path(video_in)
        )

    output:
        tuple(
            val(meta),
            path("session_${meta.id}", type: "dir")
        )

    script:
        """
        export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib"

        idtrackerai \\
            --track \\
            --output_dir . \\
            --name ${meta.id} \\
            --video_paths ${video_in} \\
            --use_bkg ${meta.bgsub} \\
            --tracking_intervals "0,${meta.video_length}" \\
            --number_of_animals 2 \\
            --intensity_ths ${meta.intensity_floor} ${meta.intensity_ceiling} \\
            --area_ths ${meta.area_floor} ${meta.area_ceiling}
        
        if [ ! -f session_${meta.id}/trajectories/with_gaps_csv/trajectories.csv ]; then
            exit 1
        fi
        """
}

process assign_ref_test {
    // determine which blob is the reference fish
    label "r_tidyverse_datatable"

    input:
        tuple(
            val(meta),
            path(session_folder)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}_traj.csv.gz")
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")

        ref <- "${meta.ref}"
        test <- "${meta.test}"
        cab_coords <- ${meta.cab_coords}

        traj_file_wo_gaps <- "${session_folder}/trajectories/without_gaps_csv/trajectories.csv"
        traj_file_with_gaps <- "${session_folder}/trajectories/with_gaps_csv/trajectories.csv"
        if (exists(traj_file_wo_gaps)) {
            message("Using trajectories without gaps")
            traj <- fread(traj_file_wo_gaps)
        } else {
            message("Using trajectories with gaps")
            traj <- fread(traj_file_with_gaps)
        }

        first_frame <- traj[1,]
        if (ref == test) {
            message("Ref is the same as test")
            idx <- 1
        } else {
            if (cab_coords == "Top") {
                idx <- which.min(first_frame[["y1"]], first_frame[["y2"]])
            } else if (cab_coords == "Bottom")
                idx <- which.max(first_frame[["y1"]], first_frame[["y2"]])
            } else if (cab_coords == "Left")
                idx <- which.min(first_frame[["x1"]], first_frame[["x2"]])
            } else if (cab_coords == "Right")
                idx <- which.max(first_frame[["x1"]], first_frame[["x2"]])
            } else {
                stop("Unrecognised ref position")
            }
        }

        if (idx == 1) {
            colnames(traj)[colnames(traj) == "x1"] <- "ref_x"
            colnames(traj)[colnames(traj) == "x2"] <- "test_x"
            colnames(traj)[colnames(traj) == "y1"] <- "ref_y"
            colnames(traj)[colnames(traj) == "y2"] <- "test_y"
        } else if (idx == 2) {
            colnames(traj)[colnames(traj) == "x1"] <- "test_x"
            colnames(traj)[colnames(traj) == "x2"] <- "ref_x"
            colnames(traj)[colnames(traj) == "y1"] <- "test_y"
            colnames(traj)[colnames(traj) == "y2"] <- "ref_y"
        } else {
            stop("Index not recognised")
        }

        traj[, fps := ${meta.fps}]
        traj[, fps := ${meta.fps}]

        fwrite(traj, "${meta.id}_traj.csv.gz")
        """
}

process visualise_identities {
    // check that ref and test are assigned correctly
    label "python_opencv_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(session_folder),
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
        font_size = 1
        font_width = 3
        colors = dict(
            ref = (0, 0, 255), # BGR red
            test = (0, 0, 0), # BGR black
        )
        font = cv2.FONT_HERSHEY_SIMPLEX
        frame = np.load("${session_folder}")
        df = pd.read_csv("${traj}", nrows = 1)
        ref_x = df.iloc[0]["ref_x"]
        ref_y = df.iloc[0]["ref_y"]
        test_x = df.iloc[0]["test_x"]
        test_y = df.iloc[0]["test_y"]

        # add points showing the test and reference fish
        cv2.circle(
            frame, tuple(ref_x, ref_y), circle_size, colors["ref"], -1
        )
        cv2.circle(
            frame, tuple(test_x, test_y), circle_size, colors["test"], -1
        )
        cv2.putText(
            frame,
            "ref",
            tuple(ref_x, ref_y),
            font,
            font_size,
            colors["ref"],
            font_width,
        )
        cv2.putText(
            frame,
            "test",
            tuple(test_x, test_y),
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


workflow TRACKING {
    take:
        split_vids
    
    main:
        track_video ( split_vids )
        //assign_ref_test ( track_video.out )
}