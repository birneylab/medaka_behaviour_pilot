#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Medaka behaviour pilot - preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

process adjust_orientation {
    // rotate 180 degrees the videos
    // camera was in the wrong orientation for the R videos
    label "ffmpeg"

    input:
        tuple(
            val(meta),
            path(video_in)
        )

    output:
        tuple(
            val(meta),
            path("${video_in.baseName}_flipped.avi")
        )

    script:
        """
        ffmpeg \\
            -i ${video_in} \\
            -q:v 1 \\
            -vf "transpose=1,transpose=1" \\
            ${video_in.baseName}_flipped.avi
        """
}

process set_split_coords {
    // visualise the splitting coordinates on a frame grab
    label "python_opencv_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(video_in)
        )

    output:
        tuple(
            val(meta),
            path("${video_in.baseName}_${meta.assay}_split_coords.png")
        )

    script:
        """
        #!/usr/bin/env python3

        import cv2 as cv
        import numpy as np

        cap = cv.VideoCapture("${video_in}")
        start = ${meta.start_frame}
        end = ${meta.end_frame}
        adj_top = ${meta.adj_top}
        adj_right = ${meta.adj_right}
        w = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
        h = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
        mid_x = round(((w - 1) / 2) + adj_right)
        mid_y = round(((h - 1) / 2) + adj_top)
        cap.set(cv.CAP_PROP_POS_FRAMES, start)
        ret, frame = cap.read()

        # Add vertical line 
        start_point = (mid_x, 0)
        end_point = (mid_x, h)
        color = (255,0,0)
        thickness = 1
        frame = cv.line(frame, start_point, end_point, color, thickness)

        # Add horizontal line
        start_point = (0, mid_y)
        end_point = (w, mid_y)
        color = (255,0,0)
        thickness = 1
        frame = cv.line(frame, start_point, end_point, color, thickness)

        # Write frame
        cv.imwrite("${video_in.baseName}_${meta.assay}_split_coords.png", frame)
        cap.release()
        """
}

process split_videos {
    // visualise the splitting coordinates on a frame grab
    label "python_opencv_numpy_pandas"

    input:
        tuple(
            val(meta),
            path(video_in)
        )

    output:
        tuple(
            val(meta),
            path("${video_in.baseName}_flipped.avi")
        )

    script:
        """
        #!/usr/bin/env python3

        import cv2 as cv
        import numpy as np

        cap = cv.VideoCapture("${video_in}")
        start = ${meta.start_frame}
        end = ${meta.end_frame}
        adj_top = ${meta.adj_top}
        adj_right = ${meta.adj_right}
        quadrant = ${meta.quadrant}
        w = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
        h = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
        mid_x = round(((w - 1) / 2) + adj_right)
        mid_y = round(((h - 1) / 2) + adj_top)
        vid_len = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
        fps = int(cap.get(cv.CAP_PROP_FPS))

        if quadrant == 'q1':
            top = 0
            bottom = mid_y
            left = mid_x
            right = w - 1
        elif quadrant == 'q2':
            top = 0
            bottom = mid_y
            left = 0
            right = mid_x
        elif quadrant == 'q3':
            top = mid_y
            bottom = h - 1
            left = 0
            right = mid_x
        elif quadrant == 'q4':
            top = mid_y
            bottom = h - 1
            left = mid_x
            right = w  - 1
        else:
            print('Invalid quadrant') 

        size = (right - left, bottom - top)

        # Define the codec and create VideoWriter object

        fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')
        out = cv.VideoWriter(OUT_FILE, fourcc, fps, size, isColor=True)

        # Capture frame-by-frame

        i = start
        while i in range(start,end):
            cap.set(cv.CAP_PROP_POS_FRAMES, i)
            # Capture frame-by-frame
            ret, frame = cap.read()
            # if frame is read correctly ret is True
            if not ret:
                print("Can't receive frame (stream end?). Exiting ...")
                break
            # Crop frame
            frame = frame[top:bottom, left:right]
            # Write frame
            out.write(frame)
            # Add to counter
            i += 1
            # Press 'esc' to close video
        #    if cv.waitKey(1) == 27:
        #        cv.destroyAllWindows()
        #        cv.waitKey(1)
        #        break

        cap.release()
        out.release()
        """
}

workflow PREPROCESSING {
    take:
        in_ch
    
    main:
        in_ch.branch {
            meta, vid ->
            L: meta.tank_side == "L"
            R: meta.tank_side == "R"
        }
        .set { in_ch_branched }
        adjust_orientation ( in_ch_branched.R )
        in_ch_branched.L
        .mix ( adjust_orientation.out )
        .set { videos_adjusted }

        videos_adjusted. combine ( ["of", "no"] )
        .map {
            meta, vid, assay ->
            [
                [
                    id: meta.id,
                    date: meta.date,
                    tank_side: meta.tank_side,
                    assay: assay,
                    fps: meta.fps,
                    // take the right variables for open field and novel object assays
                    start_frame: assay == "of" ? meta.of_start : meta.no_start,
                    end_frame: assay == "of" ? meta.of_end : meta.no_end,
                    adj_top: assay == "of" ? meta.adj_top_of : meta.adj_top_no,
                    adj_right: assay == "of" ? meta.adj_right_of : meta.adj_right_no,
                    video_length: assay == "of" ? meta.of_video_length : meta.no_video_length,
                    bgsub_q1: assay == "of" ? meta.bgsub_of_q1 : meta.bgsub_no_q1,
                    bgsub_q2: assay == "of" ? meta.bgsub_of_q2 : meta.bgsub_no_q2,
                    bgsub_q3: assay == "of" ? meta.bgsub_of_q3 : meta.bgsub_no_q3,
                    bgsub_q4: assay == "of" ? meta.bgsub_of_q4 : meta.bgsub_no_q4,
                    intensity_floor_q1: assay == "of" ? meta.intensity_floor_of_q1 : meta.intensity_floor_no_q1,
                    intensity_floor_q2: assay == "of" ? meta.intensity_floor_of_q2 : meta.intensity_floor_no_q2,
                    intensity_floor_q3: assay == "of" ? meta.intensity_floor_of_q3 : meta.intensity_floor_no_q3,
                    intensity_floor_q4: assay == "of" ? meta.intensity_floor_of_q4 : meta.intensity_floor_no_q4,
                    intensity_ceiling_q1: assay == "of" ? meta.intensity_ceiling_of_q1 : meta.intensity_ceiling_no_q1,
                    intensity_ceiling_q2: assay == "of" ? meta.intensity_ceiling_of_q2 : meta.intensity_ceiling_no_q2,
                    intensity_ceiling_q3: assay == "of" ? meta.intensity_ceiling_of_q3 : meta.intensity_ceiling_no_q3,
                    intensity_ceiling_q4: assay == "of" ? meta.intensity_ceiling_of_q4 : meta.intensity_ceiling_no_q4,
                    area_floor_q1: assay == "of" ? meta.area_floor_of_q1 : meta.area_floor_no_q1,
                    area_floor_q2: assay == "of" ? meta.area_floor_of_q2 : meta.area_floor_no_q2,
                    area_floor_q3: assay == "of" ? meta.area_floor_of_q3 : meta.area_floor_no_q3,
                    area_floor_q4: assay == "of" ? meta.area_floor_of_q4 : meta.area_floor_no_q4,
                    area_ceiling_q1: assay == "of" ? meta.area_ceiling_of_q1 : meta.area_ceiling_no_q1,
                    area_ceiling_q2: assay == "of" ? meta.area_ceiling_of_q2 : meta.area_ceiling_no_q2,
                    area_ceiling_q3: assay == "of" ? meta.area_ceiling_of_q3 : meta.area_ceiling_no_q3,
                    area_ceiling_q4: assay == "of" ? meta.area_ceiling_of_q4 : meta.area_ceiling_no_q4,
                    cab_coords_q1: assay == "of" ? meta.cab_coords_of_q1 : meta.cab_coords_no_q1,
                    cab_coords_q2: assay == "of" ? meta.cab_coords_of_q2 : meta.cab_coords_no_q2,
                    cab_coords_q3: assay == "of" ? meta.cab_coords_of_q3 : meta.cab_coords_no_q3,
                    cab_coords_q4: assay == "of" ? meta.cab_coords_of_q4 : meta.cab_coords_no_q4,
                ],
                vid
            ]
        }
        .view()
        .set { no_of_vids }
        set_split_coords ( no_of_vids )
}