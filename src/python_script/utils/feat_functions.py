"""
A
"""
import cv2
import numpy as np
from utils.vid_process_tools import color_depth_function


def draw_keypoints(video_frame, kp):
    """
    Draws keypoints onto video frame.

    Parameters
    ----------
    video_frame: ndarray
        RGB video frame input.

    kp: keypoints
        Keypoints of objects of interest.
    """
    output_frame = cv2.drawKeypoints(
        video_frame, kp, None, color=(0, 255, 0))
    return output_frame


def draw_keypoints2(video_frame, total_keypoints, depth_frame,
                    depth_camera_mtx, rot_stereo, trans_stereo):
    """
    Draws keypoints onto video frame. Also calculates X, Y, Z coordinates of
    said keypoints.

    Parameters
    ----------
    video_frame: np.ndarray
        RGB camera input.
    total_keypoints: keypoints
        List of keypoints to use.
    depth_frame: np.ndarray
        Depth camera input. Type: 16bit image.
    depth_camera_mtx: np.ndarray
        Depth camera matrix.
    rot_stereo: np.ndarray
        Stereo Rotation matrix.
    trans_stereo: np.ndarray
        Stereo Traslation matrix.

    Returns
    -------
    output_frame: np.ndarray
        RGB camera frame with drawn keypoints.
    """
    # Text configurations
    font = cv2.FONT_HERSHEY_SIMPLEX
    fontScale = 0.6
    debug_flag = True

    # Draw keypoints onto video frame.
    output_frame = cv2.drawKeypoints(video_frame, total_keypoints,
                                     None, color=(0, 0, 255))

    # Gets center of image
    (x_ct, y_ct, z_ct) = get_center_img(depth_frame, depth_camera_mtx,
                                        rot_stereo, trans_stereo)

    output_frame = cv2.circle(output_frame,
                              (int(depth_frame.shape[0]/2),
                               int(depth_frame.shape[1]/2)),
                              3,
                              (255, 255, 255),
                              3)

    # Draws 3D Coordinates on RGB frame
    cv2.putText(output_frame,
                f"({x_ct:.2f}, {y_ct:.2f}, {z_ct:.2f})",
                (int(int(depth_frame.shape[0]/2)) - 3,
                 int(int(depth_frame.shape[1]/2)) - 3), font,
                fontScale, color=(255, 255, 255), thickness=1)

    # Calculates 3D world coordinates of keypoints given known depth.
    for keypoint in total_keypoints:
        # Gets u, v coordinate of RGB image
        u_pt = keypoint.pt[0]
        v_pt = keypoint.pt[1]

        # Transforms 2D points into 3D points using depth information
        x_obj, y_obj, z_obj = transform3d(v_pt, u_pt,
                                          depth_frame, depth_camera_mtx,
                                          rot_stereo, trans_stereo)

        # Codes text  given distance
        if 1 > z_obj > 0:
            color_value = (255, 0, 0)
        elif 2 > z_obj >= 1:
            color_value = (0, 255, 0)
        else:
            color_value = (0, 0, 255)

        # Draws 3D Coordinates on RGB frame
        cv2.putText(output_frame,
                    f"({x_obj:.2f}, {y_obj:.2f}, {z_obj:.2f})",
                    (int(u_pt) - 3, int(v_pt) - 3), font,
                    fontScale, color=color_value, thickness=1)

        # Prints debugging information:
        if debug_flag:
            print("X [px]: ", u_pt, "Y [px]: ", v_pt)
            print("X [m]: ", x_obj, "Y [m]: ", y_obj, "Z [m]: ", z_obj)
    if debug_flag:
        print("Cantidad de keypoints: ", len(total_keypoints))
    return output_frame


def transform3d(u_pt, v_pt, depth, depth_camera_mtx, rot_stereo,
                trans_stereo):
    """
    Transforms (u, v) pixel coordinates into 3D coordinates given depth values.

    Parameters
    ----------
    u_pt: int
        Horizontal pixel coordinate.
    v_pt: int
        Vertical pixel coordinate
    depth_camera_mtx: np.ndarray
        Depth camera matrix.
    rot_stereo: np.ndarray
        Stereo Rotation matrix.
    trans_stereo: np.ndarray
        Stereo Traslation matrix.

    Returns
    -------
    (x_obj, y_obj, z_obj): tuple
        Tuple of 3D coordinates.
    """
    # Debugging flag
    debug_flag = True

    # Obtain (u, v, z)
    print("pixeles (u, v): ", u_pt, v_pt)
    u_pt, v_pt = int(u_pt), int(v_pt)
    z_raw = depth[u_pt][v_pt]

    # Obtain raw (x, y, z)
    z = get_distance_meters(depth[u_pt][v_pt]) / depth_camera_mtx[2][2]
    x_obj = (u_pt - depth_camera_mtx[0][2]) * z / depth_camera_mtx[0][0]
    y_obj = (v_pt - depth_camera_mtx[1][2]) * z / depth_camera_mtx[1][1]

    # Using camera matrix to proper align (x, y, z)
    (x_obj, y_obj, z_obj) = (rot_stereo.dot((x_obj, y_obj, z))
                             + trans_stereo)

    # Prints debugging information
    if debug_flag:
        print("z_distance[px]: ", depth[u_pt][v_pt])
        print("z_distance[m]: ", get_distance_meters(depth[u_pt][v_pt]))
        print("output z[m]: ", z_obj)

    # Output mode:
    output_select = "z_in_meters"
    output_dict_mode = {
        "normal": (x_obj, y_obj, z_obj),
        "z_in_meters": (x_obj, y_obj, z),
        "raw_z_value": (x_obj, y_obj, z_raw),
        "k_coord_2_std_coord": (-x_obj, z_obj, y_obj),
    }
    return output_dict_mode[output_select]


def get_center_img(depth, depth_camera_mtx, rot_stereo,
                   trans_stereo):
    """
    Registers center of depth image onto 3D coordinates

    Parameters
    ----------
    depth_camera_mtx: np.ndarray
        Depth camera matrix.
    rot_stereo: np.ndarray
        Stereo Rotation matrix.
    trans_stereo: np.ndarray
        Stereo Traslation matrix.

    Returns
    -------
    (x_obj, y_obj, z_obj): tuple
        Tuple of 3D coordinates.
    """
    # Obtain image center representation in 3D
    depth_shape = depth.shape
    u_ct, v_ct = int(depth_shape[0]/2), int(depth_shape[1]/2)
    print(u_ct, v_ct)
    z_center = get_distance_meters(
        depth[u_ct][v_ct]) / depth_camera_mtx[2][2]
    x_ct = (u_ct - depth_camera_mtx[0][2]) * z_center / depth_camera_mtx[0][0]
    y_ct = (v_ct - depth_camera_mtx[1][2]) * z_center / depth_camera_mtx[1][1]

    # Using camera matrix to proper align (x, y, z)
    (x_obj, y_obj, z_obj) = (rot_stereo.dot((x_ct, y_ct, z_center))
                             + trans_stereo)

    return (x_obj, y_obj, z_obj)


def get_distance_meters(measurement_obj):
    """
    Get distance in meters of depth measurement.

    Parameters
    ----------
    measurement_obj: float
        Measurement of depth. Must be lower than 2047.
    """
    linear_approx = (measurement_obj * -0.0030711016) + 3.3309495161
    output_measurement = 1.0 / linear_approx
    return output_measurement


def mark_point_interest(x_pt, y_pt, depth_frame, binary_mask):
    """
    Draws certain point on screen and shows its distance.
    """
    # Get depth measurement of point of interest
    depth_measurement_meters = get_distance_meters(depth_frame[x_pt, y_pt])
    multipl = np.multiply(binary_mask.astype(np.uint8)/255, depth_frame)

    heatmap = color_depth_function(multipl)
    binary_mask_0_1 = binary_mask.astype(np.uint8)/255
    binary_mask_3D = np.dstack((
        binary_mask_0_1, binary_mask_0_1, binary_mask_0_1))
    heatmap = np.multiply(binary_mask_3D, heatmap)

    # Get measurements of point of interest
    str_measurement_m = "Measurement [m]: " + str(depth_measurement_meters)
    str_measurement_px = "Measurement [pixel value]: " + str(
        depth_frame[x_pt, y_pt])

    # Draw point and show text on screen
    cv2.circle(depth_frame, (x_pt, y_pt), 2, (255, 255, 255), 2)
    cv2.circle(heatmap, (x_pt, y_pt), 2, (255, 255, 255), 2)
    font = cv2.FONT_HERSHEY_SIMPLEX
    fontScale = 0.5
    cv2.putText(heatmap, str_measurement_m, (x_pt - 3, y_pt - 3), font,
                fontScale, color=(255, 0, 0), thickness=2)
    cv2.putText(heatmap, str_measurement_px, (x_pt - 3, y_pt - 15),
                font, fontScale, color=(255, 0, 0), thickness=2)
    return heatmap, depth_frame


def convert(img, target_type_min, target_type_max, target_type):
    imin = img.min()
    imax = img.max()

    a = (target_type_max - target_type_min) / (imax - imin)
    b = target_type_max - a * imax
    new_img = (a * img + b).astype(target_type)
    return new_img


def registerDepth(intrinsics_depth, intrinsics_color, dist_coef_color,
                  extrinsics, depth, shape):
    assert dist_coef_color is None

    width, height = shape
    out = np.zeros((height, width))/0
    y, x = np.meshgrid(np.arange(depth.shape[0]),
                       np.arange(depth.shape[1]), indexing='ij')
    x = x.reshape(1, -1)
    y = y.reshape(1, -1)
    z = depth.reshape(1, -1)
    x = (x - intrinsics_depth[0, 2])/intrinsics_depth[0, 0]
    y = (y-intrinsics_depth[1, 2])/intrinsics_depth[1, 1]
    pts = np.vstack((x*z, y*z, z))
    pts = extrinsics[:3, :3]@pts+extrinsics[:3, 3:]
    pts = intrinsics_color@pts
    px = np.round(pts[0, :]/pts[2, :])
    py = np.round(pts[1, :]/pts[2, :])
    mask = (px >= 0) * (py >= 0) * (px < width) * (py < height)
    out[py[mask].astype(int), px[mask].astype(int)] = pts[2, mask]
    return out


def registerDepth2(intrinsics_depth, intrinsics_color, dist_coef_color,
                   rot_mtx, trans_vect, depth, shape):
    assert dist_coef_color is not None

    width, height = shape
    out = np.zeros((height, width))/0
    y, x = np.meshgrid(np.arange(depth.shape[0]),
                       np.arange(depth.shape[1]), indexing='ij')
    x = x.reshape(1, -1)
    y = y.reshape(1, -1)
    z = depth.reshape(1, -1)
    print(intrinsics_depth[0])
    x = (x - intrinsics_depth[0][2])/intrinsics_depth[0][0]
    y = (y-intrinsics_depth[1][2])/intrinsics_depth[1][1]
    pts = np.vstack((x*z, y*z, z))
    pts = rot_mtx@pts+trans_vect
    pts = intrinsics_color@pts
    px = np.round(pts[0, :]/pts[2, :])
    py = np.round(pts[1, :]/pts[2, :])
    mask = (px >= 0) * (py >= 0) * (px < width) * (py < height)
    out[py[mask].astype(int), px[mask].astype(int)] = pts[2, mask]
    return out
