"""
A
"""
import numpy as np
import cv2
import freenect
import skimage.exposure


def get_video(corner_detector, mtx_obj, dist_obj, frame_size, ORB=True):
    """
    Get video frame and apply ORB/FAST corner detector.
    """
    # Get video frame
    video_obj = freenect.sync_get_video()[0]
    video_obj = cv2.cvtColor(video_obj, cv2.COLOR_BGR2RGB)

    # Undistorts image
    undistort = False
    if undistort:
        video_obj = undistort_img(video_obj, mtx_obj, dist_obj, frame_size)

    if ORB:
        # find the keypoints with ORB
        kp = corner_detector.detect(video_obj, None)
        # compute the descriptors with ORB
        kp, _ = corner_detector.compute(video_obj, kp)

    else:
        # Apply FAST Corner detector
        kp = corner_detector.detect(video_obj, None)

    return video_obj, kp


def get_infrared():
    """
    Retrieves infrared information
    """
    array, _ = freenect.sync_get_video(0, freenect.VIDEO_IR_10BIT)
    return array


def get_depth(mtx_obj, dist_obj, frame_size):
    """
    get depth frame.
    """
    # Calculate depth frame. Returns 16-bit image.
    depth_obj = freenect.sync_get_depth()[0]
    # Undistorts image
    undistort = False
    if undistort:
        depth_obj = undistort_img(depth_obj, mtx_obj, dist_obj, frame_size)

    return depth_obj


def undistort_img(img_obj, mtx_obj, dist_obj, frame_size):
    """
    Undistorts image
    """
    new_mtx_obj, roi = cv2.getOptimalNewCameraMatrix(
        mtx_obj, dist_obj, frame_size, 1, frame_size)
    output = cv2.undistort(img_obj, mtx_obj, dist_obj, None, new_mtx_obj)
    x, y, w, h = roi
    output = output[y:y+h, x:x+w]
    return output


def pretty_depth(depth):
    np.clip(depth, 0, 2**10-1, depth)
    depth >>= 2
    depth = depth.astype(np.uint8)
    return depth


def color_depth_function(img_obj):
    """
    a
    """
    # stretch to full dynamic range
    stretch = skimage.exposure.rescale_intensity(
        img_obj, in_range='image', out_range=(0, 255)).astype(np.uint8)

    # convert to 3 channels
    stretch = cv2.merge([stretch, stretch, stretch])

    # define colors
    color1 = (0, 0, 255)     # red
    color2 = (0, 165, 255)   # orange
    color3 = (0, 255, 255)   # yellow
    color4 = (255, 255, 0)   # cyan
    color5 = (255, 0, 0)     # blue
    color6 = (128, 64, 64)   # violet
    colorArr = np.array(
        [[color1, color2, color3, color4, color5, color6]], dtype=np.uint8)

    # resize lut to 256 (or more) values
    lut = cv2.resize(colorArr, (256, 1), interpolation=cv2.INTER_LINEAR)

    # apply lut
    result = cv2.LUT(stretch, lut)

    # create gradient image
    grad = np.linspace(10, 255, 512, dtype=np.uint8)
    grad = np.tile(grad, (20, 1))
    grad = cv2.merge([grad, grad, grad])
    _ = cv2.LUT(grad, lut)
    return result


def depth_to_heatmap(depth_frame):
    """
    Converts depth frame to heatmap, to visualize it better.

    Parameters
    ----------
    depth_frame: np.ndarray

    """
    depth_frame = cv2.normalize(
        depth_frame, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U
        )
    depth_frame = cv2.bitwise_not(depth_frame)
    depth_frame = cv2.applyColorMap(depth_frame, cv2.COLORMAP_JET)
    return depth_frame


def smooth(infrared_frame):
    infrared_frame = cv2.medianBlur(infrared_frame, 19)
    infrared_frame = cv2.bilateralFilter(infrared_frame, 9, 75, 75)
    infrared_frame = cv2.GaussianBlur(infrared_frame, (7, 7), 0)
    return infrared_frame


def fill_misdetections_depth(depth_frame_uint16):
    """
    Fill misdetections in depth image.

    Parameters
    ----------
    depth_frame_uint16: np.uint16
        Raw depth image

    Returns
    -------
    depth_impaint: np.uint16
        Impainted depth image, with patched holes.
    """
    binary_mask = threshold_mask(depth_frame_uint16)
    depth_impaint = cv2.inpaint(depth_frame_uint16,
                                binary_mask.astype(np.uint8), 3,
                                cv2.INPAINT_TELEA)
    return depth_impaint


def threshold_mask(img_obj):
    """
    a
    """
    _, mask = cv2.threshold(img_obj, 2046, 2047, cv2.THRESH_BINARY)
    return mask
