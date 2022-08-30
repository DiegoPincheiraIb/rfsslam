"""
realt_class.py
    Class that work with video.
"""
import numpy as np
import sys
import freenect
import cv2
import copy
import os
import time
import yaml
from datetime import datetime
from utils.misc_tools import create_folder


class VideoDisplay():
    """
    Class that operates by real time video.
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
            "depth": None,
            "ir": None,
        }
        self.help_dict = {
            "h": "Shows help commands.",
            "c": "Saves onto disk current RGB, Depth and IR frame.",
            "q": "Quits program.",
            "v": "Records video of RGB, Depth and IR input."
        }

        self.flags = {
            "is recording": False,
            "was initial tstmp folder created": True
        }

        self.video_objects = {
            "rgb": None,
            "depth": None,
            "ir": None,
        }

        self.paths = {
            "main path": "data/rgbd/",
            "initial timestamp": copy.deepcopy(
                datetime.now().strftime("%d_%b_%Y_%H_%M_%S")),
            "current timestamp": None,
        }
        # Imports main configuration file
        self.main_yaml = main_yaml

    def run_algorithms(self):
        """
        Run main algorithms.
        """
        # Begin main loop
        while True:
            # Define key press variable
            key_press = cv2.waitKey(1)

            self.show_help(key_press)

            # Main algorithms
            self.get_video()
            self.get_depth()
            self.get_ir()
            for video_type in ["rgb", "depth", "ir"]:
                print(
                    video_type + " dtype: ", self.img_cache[video_type].dtype)

            # Begin capture video
            self.begin_capture_video(key_press)
            if self.flags["is recording"]:
                self.capture_video()
            self.stop_capture_video(key_press)

            # Shows frames
            self.show_frames()

            # Saves frames onto disk
            self.save_frames(key_press)

            # Closes program
            self.close_program(key_press)

    def get_video(self):
        """
        Gets video frame onto image cache.
        """
        # Get video frame
        video_obj = freenect.sync_get_video()[0]
        video_obj = cv2.cvtColor(video_obj, cv2.COLOR_BGR2RGB)
        self.img_cache["rgb"] = video_obj

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

    def save_frames(self, key_press):
        """
        Save frames onto disk.
        """
        if key_press == ord("c"):
            timestamp_obj = datetime.now().strftime(
                "%d_%b_%Y_%H_%M_%S")
            # Create folders (omitted if folder was created)
            path_tstmp_folder = create_folder(
                    self.paths["main path"], self.paths["initial timestamp"]
            )
            path_imgs_folder = create_folder(
                    path_tstmp_folder, "images"
            )
            for str_type in ["rgb", "depth", "ir"]:
                _ = create_folder(
                        path_imgs_folder, str_type
                )
            for str_type in list(self.img_cache.keys()):
                # Saves frame with corresponding type.
                # E.g. frame_rgb_23847234.npy corresponds to a rgb image.
                np.save(
                    path_imgs_folder + "/"
                    + str_type + "/"
                    + "frame_" + str_type + "_" + str(timestamp_obj) + ".npy",
                    self.img_cache[str_type])

                # Rewrites yaml with new information
                # self.rewrite_yaml(timestamp_obj)

    def rewrite_yaml(self, timestamp_obj):
        """
        Writes information of saved frames onto yaml.
        """
        str_tstamp = str(timestamp_obj)

        # Checks if timestamp is not on list
        if str_tstamp not in self.main_yaml["frames"]["list_timestamps"]:

            # Add timestamp onto main yaml
            self.main_yaml["frames"]["list_timestamps"].append(
                str(timestamp_obj))

            # Rewrite yaml onto disk
            with open(
                    "src/python_script/config/cfg_parameters.yaml",
                    "w", encoding="utf-8") as yaml_file:
                yaml_file.write(yaml.dump(self.main_yaml,
                                          default_flow_style=False))

    def begin_capture_video(self, key_press):
        """
        Signals program to begin capturing information from input sources
        (rgb, depth and IR).
        Saves them onto disk.
        """
        if key_press == ord("v") and not self.flags["is recording"]:
            total_countdown = 5
            for idx in range(total_countdown):
                print('\007')
                print("Beginning recording in ", total_countdown - idx)
                time.sleep(1)
            print("Recording has started!")
            self.flags["is recording"] = True

            # Registers time of initialization of recording in "timestamp_obj"
            current_time = datetime.now()
            self.paths["current timestamp"] = current_time.strftime(
                "%d_%b_%Y_%H_%M_%S")

            # Initializes output video folders
            self.paths["current_data_dir"] = (
                self.paths["main path"]
                + self.paths["current timestamp"] + "/")
            self.paths["current_img_dir"] = (
                self.paths["current_data_dir"] + "images/")
            # os.mkdir(self.paths["main path"])
            os.mkdir(self.paths["current_data_dir"])
            os.mkdir(self.paths["current_img_dir"])
            for video_type in ["rgb", "depth", "ir"]:
                os.mkdir(self.paths["current_img_dir"] + video_type)

    def capture_video(self):
        """
        Saves input video stream.
        """
        current_frame = datetime.now()
        str_current_frame = current_frame.strftime(
            "%d_%b_%Y_%H_%M_%S")
        for video_type in ["rgb", "depth", "ir"]:
            # Saves frame with corresponding type.
            # E.g. frame_rgb_23847234.npy corresponds to a rgb image.
            np.save(
                self.paths["current_img_dir"] + video_type + "/"
                + "frame_" + video_type + "_"
                + str(str_current_frame) + ".npy",
                self.img_cache[video_type])

            # Rewrites yaml with new information
            # self.rewrite_yaml(str_current_frame)

    def stop_capture_video(self, key_press):
        """
        Stops capturing input video stream.
        """
        if key_press == ord("s") and self.flags["is recording"]:
            print("Recording has stopped!")
            self.flags["is recording"] = False

    def show_help(self, key_press):
        """
        Shows help on screen.
        """
        if key_press == ord("h"):
            for key_obj in self.help_dict:
                print(key_obj, ":", self.help_dict[key_obj])

    def show_frames(self):
        """
        Shows frames saved on image cache.
        """
        for frame_type in ["rgb", "depth", "ir"]:
            cv2.imshow(frame_type, self.img_cache[frame_type])

    def close_program(self, key_press):
        """
        Closes program if key is pressed.
        """
        # Press "q" to quit process
        if key_press == ord("q"):
            cv2.destroyAllWindows()
            sys.exit(1)
