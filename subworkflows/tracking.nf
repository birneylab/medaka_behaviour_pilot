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
        export LD_LIBRARY_PATH="\$CONDA_PREFIX/lib"

        idtrackerai \\
            --video_paths ${video_in} \\
            --output_dir . \\
            --use_bkg ${meta.bgsub} \\
            --tracking_intervals "0,${meta.video_length}" \\
            --number_of_animals 2 \\
            --intensity_ths ${meta.intensity_floor} ${meta.intensity_ceiling} \\
            --area_ths ${meta.area_floor} ${meta.area_ceiling} \\
            --track
        
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

workflow TRACKING {
    take:
        split_vids
    
    main:
        track_video ( split_vids )
}