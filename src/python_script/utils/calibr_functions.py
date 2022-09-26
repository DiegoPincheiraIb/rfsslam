"""
A
"""
import cv2
import numpy as np


def retrieve_corners(img_gray, checkboard_dim=(8, 6)):
    """
    Calculates object points and image points.
    """
    ret_obj, corners_obj = cv2.findChessboardCorners(
        img_gray, checkboard_dim, None)
    return ret_obj, corners_obj


def add_objp_imgp(img_obj, img_gray, ret_obj, corners_obj, objpoints,
                  imgpoints, criteria, objp, key_press, checkboard_dim=(8, 6)):
    """
    Add image points and object points
    """
    if key_press == ord("c"):
        # If found, add object points, image points (after refining them)
        corners_obj = cv2.cornerSubPix(img_gray, corners_obj, (11, 11),
                                       (-1, -1), criteria)
        imgpoints.append(corners_obj)
    # Draw and display the corners
    cv2.drawChessboardCorners(img_obj, checkboard_dim, corners_obj, ret_obj)
    return img_obj, objpoints, imgpoints


def generate_calibration_parameters():
    """
    Readies calibration parameters
    """
    chessboard_size = (8, 6)

    # termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER,
                30, 0.001)

    # prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
    objp = np.zeros((chessboard_size[0] * chessboard_size[1], 3), np.float32)
    objp[:, :2] = np.mgrid[0:chessboard_size[0],
                           0:chessboard_size[1]].T.reshape(-1, 2)
    objpoints = []  # 3d point in real world space
    imgpointsL = []  # 2d points in image plane.
    imgpointsR = []  # 2d points in image plane.
    return objpoints, imgpointsL, imgpointsR, objp, criteria, chessboard_size


def calibrate_one_camera(
        objpoints, imgpoints_obj, frame_size, img_obj):
    """
    a
    """
    ret_obj, camera_matrix, dist_obj, rvecs, tvecs = cv2.calibrateCamera(
        objpoints, imgpoints_obj, frame_size, None, None)
    height_obj, width_obj, _ = img_obj.shape
    new_camera_matrix, _ = cv2.getOptimalNewCameraMatrix(
        camera_matrix, dist_obj, (width_obj, height_obj), 1,
        (width_obj, height_obj))
    return new_camera_matrix, dist_obj, rvecs, tvecs, ret_obj


def calibrate_stereo(
        objpoints, imgpointsL, imgpointsR, new_cam_mtx_L, distL,
        new_cam_mtx_R, distR, grayL):
    """
    a
    """
    flags = 0
    flags |= cv2.CALIB_FIX_INTRINSIC
    # Here we fix the intrinsic camara matrixes so that only Rot, Trns,
    # Emat and Fmat are calculated.
    # Hence intrinsic parameters are the same
    criteria_stereo = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER,
                       30, 0.001)
    (retStereo, new_cam_mtx_L,
     distL, new_cam_mtx_R, distR,
     rot, trans, essential_matrix,
     fundamentalMatrix) = cv2.stereoCalibrate(
         objpoints, imgpointsL, imgpointsR, new_cam_mtx_L, distL,
         new_cam_mtx_R, distR, grayL.shape[::-1], criteria_stereo, flags)

    # Return parameters
    return (retStereo, new_cam_mtx_L,
            distL, new_cam_mtx_R, distR,
            rot, trans, essential_matrix,
            fundamentalMatrix)
