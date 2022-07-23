"""
realt_class.py
    Class that work with video.
"""
import numpy as np
import sys
import freenect
import cv2
import copy


class CalibrDisplay():
    """
    Class that calibrates camera.
    """
    def __init__(self, main_yaml):
        """
            Initialization method.

            Parameters
            ----------

            Returns
            -------
            None.
        """
        # Initializes image cache.
        self.img_cache = {
            "rgb": None,
            "ir": None,
        }

        # Imports main configuration file
        self.main_yaml = main_yaml

        # Initializes calibration parameters
        self.chessbd_size = self.main_yaml["calibr"]["chessboard_size"]

        # termination criteria
        self.criteria = (cv2.TERM_CRITERIA_EPS
                         + cv2.TERM_CRITERIA_MAX_ITER,
                         30, 0.001)
        self.w_size = (self.main_yaml["calibr"]["w_size"]["x"],
                       self.main_yaml["calibr"]["w_size"]["y"])
        self.zero_zone = (self.main_yaml["calibr"]["zero_zone"]["x"],
                          self.main_yaml["calibr"]["zero_zone"]["y"])

        # prepare object points, like (0,0,0), (1,0,0), (2,0,0) ...,(6,5,0)
        self.objp = np.zeros((
            self.chessbd_size["x"] * self.chessbd_size["y"], 3),
            np.float32)
        self.objp[:, :2] = np.mgrid[0:self.chessbd_size["x"],
                                    0:self.chessbd_size["y"]
                                    ].T.reshape(-1, 2)
        self.objpoints = []  # 3d point in real world space
        self.imgpointsL = []  # 2d points in image plane.
        self.imgpointsR = []  # 2d points in image plane.

    def run_algorithms(self):
        """
        Run main algorithms
        """
        # Begin main loop
        print((self.chessbd_size["x"], self.chessbd_size["y"]))
        while True:
            # Define key press variable
            key_press = cv2.waitKey(1)
            # Main algorithms
            self.get_video()
            self.get_grayscale_img()
            self.get_ir()

            self.img_cache["ir"] = self.pretty_depth(self.img_cache["ir"])

            # Get corners of RGB image
            retL, cornersL = cv2.findChessboardCorners(
                self.img_cache["gray"],
                (self.chessbd_size["x"], self.chessbd_size["y"]), None)
            # Get corners of IR Image
            retR, cornersR = cv2.findChessboardCorners(
                self.img_cache["ir"],
                (self.chessbd_size["x"], self.chessbd_size["y"]), None)

            # If found checkerboard in both images:
            if retL:
                # RGB Image
                cornersL = self.generate_subcorner_pix("gray", cornersL)
                self.show_checkerbd("rgb", cornersL, retL)
            if retR:
                # RGB Image
                cornersR = self.generate_subcorner_pix("ir", cornersR)
                self.show_checkerbd("ir", cornersR, retR)
            # Shows frames
            self.show_frames()

            # Closes program
            self.close_program(key_press)

    def pretty_depth(self, depth):
        np.clip(depth, 0, 2**10-1, depth)
        depth >>= 2
        depth = depth.astype(np.uint8)
        return depth

    def get_video(self):
        """
        Gets video frame onto image cache.
        """
        # Get video frame
        video_obj = freenect.sync_get_video()[0]
        video_obj = cv2.cvtColor(video_obj, cv2.COLOR_BGR2RGB)
        self.img_cache["rgb"] = video_obj

    def get_grayscale_img(self):
        """
        Gets grayscale image from rgb
        """
        img_copy = copy.deepcopy(self.img_cache["rgb"])
        self.img_cache["gray"] = cv2.cvtColor(img_copy, cv2.COLOR_BGR2GRAY)

    def get_depth(self):
        """
        Gets depth frame.
        """
        # Calculate depth frame. Returns 16-bit image.
        depth_obj = freenect.sync_get_depth()[0]
        self.img_cache["depth"] = depth_obj

    def get_ir(self):
        """
        Gets infrared frame.
        """
        array, _ = freenect.sync_get_video(0, freenect.VIDEO_IR_10BIT)
        self.img_cache["ir"] = array

    def generate_subcorner_pix(self, str_mode, corners_obj):
        """
        A
        """
        corners_filtered = cv2.cornerSubPix(
            self.img_cache[str_mode], corners_obj,
            self.w_size,
            self.zero_zone,
            self.criteria)
        return corners_filtered

    def show_checkerbd(self, str_mode, corners, ret_obj):
        """
        Shows checkerboard grid on RGB or Infrared image.
        """
        cv2.drawChessboardCorners(
            self.img_cache[str_mode],
            (self.chessbd_size["x"],
             self.chessbd_size["y"]),
            corners, ret_obj)

    def capture_checkerbd(self):
        """
        Saves checkerboard info onto RAM.
        """

    def show_frames(self):
        """
        Shows frames saved on image cache.
        """
        for frame_type in ["rgb", "ir", "gray"]:
            cv2.imshow(frame_type, self.img_cache[frame_type])

    def close_program(self, key_press):
        """
        Closes program if key is pressed.
        """
        # Press "q" to quit process
        if key_press == ord("q"):
            cv2.destroyAllWindows()
            sys.exit(1)
