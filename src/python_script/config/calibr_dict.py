"""
calibr_dict.py
    Dictionary of calibration parameters.
"""

import cv2

calibr_process_dict = {
    "x_sqrs": 7,
    "y_sqrs": 6,
    "n_dims": 3,
    "criteria": (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER,
                 30, 0.001),
    "w_size": (11, 11),
    "zero_zone": (-1, -1),

}
