"""
frames_class.py
    Class that work with frames.
"""
import os
import shutil
from shutil import SameFileError
import sys
import copy
import time
import numpy as np
import cv2
import json
import pandas as pd
from datetime import datetime
from utils.ssc import ssc
from utils.misc_tools import (
    merge_R_t_one_mtx,
    load_frames, change_img
)
from utils.feat_functions import (
    draw_keypoints,
    get_distance_meters,
)


class FrameDisplay():
    """
    Class that operates frame by frame.

    Booleans
    --------
    are_3d_points_loaded
        Flag that indicates that 3d points were loaded.
    are_kp_loaded
        Flag that indicates that keypoints were loaded.
    are_keypts_shown
        Flag that indicates that keypoints are shown on screen.
    were_3dpts_calc
        Flag that indicates that 3d points were calculated, not loaded.
    is_depth_processed
        Flag to process depth image for visualization.
    """
    def __init__(self, main_yaml, dict_calibr):
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
            "nextrgb": None,
            "depth": None,
            "ir": None,
            "3dpts": None,
            "kps": None,
            "mask": None
        }

        self.flags = {
            "is recording": False,
        }

        # Imports main configuration file
        self.main_yaml = main_yaml

        # Initializes dictionary containing information of all frames.
        self.frames_dict = {}

        # Loads list of timestamps
        self.get_list_frames_from_folder()

        # Initializes current id and current timestamp
        self.current_id = 0
        self.timestamp = self.main_yaml["frames"]["list_timestamps"][0]

        self.new_timestamp = self.main_yaml["new_timestamp"]

        # Initializes amount of processed images
        self.amount_img_calc = 0

        self.list_3dpts = []

        # Initializes pandas dataframe to store pose estimation
        self.df_pose = pd.DataFrame(
            [[float(0), float(0), float(0.1),
             float(0), float(0), float(0), float(0)]],
            index=range(len(self.main_yaml["frames"]["list_timestamps"])),
            columns=["dx", "dy", "dz", "dqx", "dqy", "dqz", "dqw"])
        self.df_pose.loc[0] = [float(0) for _ in range(6)].append(float(1))

        # Set to calculate 3d points
        self.are_3d_points_loaded = False
        self.are_kp_loaded = False
        self.are_keypts_shown = False
        self.were_3dpts_calc = False

        # Flag to process depth image for visualization.
        self.is_depth_processed = False

        # Intrinsic parameters:
        # RGB
        self.rgb_camera_dict = dict_calibr["RGB"]
        self.rgb_cam_mtx = np.array(self.rgb_camera_dict["mtx"])
        self.rgb_dist = np.array(self.rgb_camera_dict["dist"]).flatten()

        # Depth
        self.depth_camera_dict = dict_calibr["Depth"]
        self.depth_cam_mtx = np.array(self.depth_camera_dict["mtx"])
        self.depth_dist = np.array(self.depth_camera_dict["dist"]).flatten()

        # Extrinsic parameters:
        # Stereo
        self.stereo_camera_dict = dict_calibr["Stereo"]
        self.trans_stereo = np.array(
            self.stereo_camera_dict["trans"]).flatten()
        self.rot_stereo = np.array(self.stereo_camera_dict["rot"])
        self.merged_r_t = merge_R_t_one_mtx(self.rot_stereo, self.trans_stereo)

        # Initializes ORB detector
        self.orb_detector = cv2.ORB_create(
            nfeatures=main_yaml["orb_features"]["n_feats"])

        # Fixes frame size
        self.frame_size = (main_yaml["camera"]["width"],
                           main_yaml["camera"]["height"])

        # Initializes 3d coordinates dictionary
        self.coord_values = {}

    def load_first_frame(self):
        """
        Loads first frame onto img_cache.
        Parameters
        ----------

        Returns
        ------
        """
        # Fills dictionary with array of frames.
        path_obj = ("data/rgbd/" + self.main_yaml["frames"]["chosen_id"]
                    + "/images/")
        self.frames_dict = load_frames(
            path_obj,
            self.main_yaml["frames"]["list_timestamps"])
        print("Largo de frames: ", len(self.frames_dict))
        # Loads onto current image cache the first frame.
        for frame_type in ["rgb", "depth", "ir"]:
            self.img_cache[frame_type] = (
                self.frames_dict[
                    self.main_yaml["frames"]["list_timestamps"][
                        self.current_id]
                    ][frame_type]
            )
            self.img_cache["nextrgb"] = (
                self.frames_dict[
                    self.main_yaml["frames"]["list_timestamps"][
                        self.current_id + 1]
                    ]["rgb"]
            )
        self.threshold_mask()
        print("Loading 3d points. Please wait...")
        # self.load_3dpts_mtx()
        print("3d points loaded from file.")

    def run_algorithms(self):
        """
        Run main algorithms.

        Parameters
        ----------

        Returns
        ------
        """
        # Begin main loop
        while True:
            # Define key press variable
            key_press = cv2.waitKey(1)

            # Changes current image
            self.change_img(key_press)

            orb_comparison_mode = False
            if orb_comparison_mode:
                cv2.imshow("aa", self.img_cache["nextrgb"])

                (queryKeypoints,
                    queryDescriptors) = self.orb_detector.detectAndCompute(
                        self.img_cache["rgb"], None)
                (trainKeypoints,
                    trainDescriptors) = self.orb_detector.detectAndCompute(
                        self.img_cache["nextrgb"], None)

                matcher = cv2.BFMatcher()
                matches = matcher.match(queryDescriptors, trainDescriptors)
                final_img = cv2.drawMatches(
                    self.img_cache["rgb"], queryKeypoints,
                    self.img_cache["nextrgb"], trainKeypoints,
                    matches[:20], None)

                final_img = cv2.resize(final_img, (1000, 650))

                # Show the final image
                cv2.imshow("Matches", final_img)

            mode_obj = "visualization"  # "3d points" # "visualization"

            if mode_obj == "3d points":

                # Gets 3d points of current frame and saves them onto RAM.
                self.get_3dpts()
                self.check_3dpts()
                # Gets keypoints of current frame and saves them onto RAM.
                self.get_keypts()
                self.check_keypts()

                # Check if keypoints are shown on video frame
                if not self.are_keypts_shown:
                    self.get_u_v_from_3dpts()
                    self.match_keypts_2dpts()
                    self.are_keypts_shown = True

                # Depth 16 bit to 8 bit
                if not self.is_depth_processed:
                    # self.fill_misdetect_depth()
                    self.invert_depth_image()
                    self.convert_depth_16b_2_8bit()
                    # self.fill_misdetect_depth()
                    self.is_depth_processed = True

                # Takes capture frame
                self.take_capture(key_press)

                # Generates 3D Array
                self.generate_3D_array(key_press)

                self.begin_measurement_process(key_press)

                if self.flags["is recording"]:
                    # self.save_3dpoints_mtx()
                    self.measurement_process()

            # input_pose_estimation
            self.input_pose_estimation(key_press)

            self.copy_current_file_to_new_timestamp(key_press)

            # Shows frames
            self.show_frames()

            # Closes program
            self.close_program(key_press)

    def change_img(self, key_press):
        """
        Changes current image.

        Parameters
        ----------

        Returns
        ------
        """
        # Rotates over images on folder.
        # Activates if there are more than one image on folder
        # And n or b keys were pressed.
        if (key_press == ord("n") or key_press == ord("b") and (
                len(self.frames_dict) > 1)):

            # Updates current id, current timestamp, and image cache.
            (self.img_cache,
                self.current_id,
                self.timestamp) = change_img(
                self.main_yaml["key_press"][str(key_press)],
                self.img_cache,
                self.current_id,
                self.main_yaml["frames"]["list_timestamps"],
                self.frames_dict,
                self.timestamp)
            self.threshold_mask()
            # Resets flags so that 3d points and keypoints can be loaded.
            self.are_3d_points_loaded = False
            self.are_kp_loaded = False
            self.are_keypts_shown = False
            self.is_depth_processed = False
            self.were_3dpts_calc = False

    def get_3dpts(self):
        """
        Gets 3d points.
        """
        # Checks if 3D points are calculated for current frame
        # If there is not 3d points associated with that frame,
        # calculates them and adds them onto img_cache.
        if "3dpts" not in list(self.frames_dict[self.timestamp].keys()):
            # Initializes 3d points matrix.
            mtx_3dpt = np.zeros((480, 640, 3))

            # Test depth
            depth_3d_v2 = np.zeros((480, 640, 3))
            depth_3d_v3 = np.zeros((480, 640, 3))

            # Calculates 3d points
            depth_3d = cv2.rgbd.depthTo3d(self.img_cache["depth"],
                                          self.depth_cam_mtx)

            # Applies extrinsic parameters to 3d points and
            # saves them onto dictionary of frames.
            for u in range(depth_3d.shape[0]):
                for v in range(depth_3d.shape[1]):
                    mtx_3dpt[u, v] = (self.rot_stereo.dot(depth_3d[u, v])
                                      + self.trans_stereo)
                    statement_obj = False
                    if statement_obj:
                        if u % 100 == 0 and v % 100 == 0:
                            # print("depth3d [u,v]", depth_3d[u, v])
                            # print("3d points: ", mtx_3dpt[u, v])
                            print("points: ", u, v)
                            print("generated by function: ", depth_3d[u, v])
                            print("rotated by function: ", mtx_3dpt[u, v])
                            print("generated by hand: ", depth_3d_v2[u, v])
                            print("rotated by function: ", depth_3d_v3[u, v])
                            print("--------")
            self.frames_dict[self.timestamp]["3dpts"] = mtx_3dpt
            print("Saved 3D points onto RAM.")

    def get_3dpts_v2(self):
        """
        Calculates 3d points from depth image.
        """
        if "3dpts" not in list(self.frames_dict[self.timestamp].keys()):
            # Converts depth from raw to meters.
            # self.convert_depth()

            # focal lengths
            x_focal = self.depth_cam_mtx[0][0]
            y_focal = self.depth_cam_mtx[1][1]

            depth_3d = copy.deepcopy(self.img_cache["depth"])
            mtx_3dpt = np.zeros((480, 640, 3))
            for u in range(depth_3d.shape[0]):
                for v in range(depth_3d.shape[1]):
                    zw = x_focal / np.sqrt(u**2 + v**2 + x_focal**2)
                    xw = zw * u / x_focal
                    yw = zw * v / y_focal
                    mtx_3dpt[u, v] = np.array([xw, yw, zw])
                    mtx_3dpt[u, v] = (
                        self.rot_stereo.dot(np.array([xw, yw, zw]))
                        + self.trans_stereo)
                    statement_obj = True
                    if statement_obj:
                        if u % 100 == 0 and v % 100 == 0:
                            print("depth3d [u,v]", depth_3d[u, v])
                            print("3d points: ", mtx_3dpt[u, v])
                            print("--------")
            self.frames_dict[self.timestamp]["3dpts"] = mtx_3dpt
            print("Saved 3D points onto RAM.")

    def get_3dpts_v3(self):
        """
        Calculates 3d points given
        """

    def check_3dpts(self):
        """
        Checks if 3d points are loaded onto image cache.
        """
        # Loads 3d points onto img cache
        if not self.are_3d_points_loaded:
            # Retrieves them from dictionary
            self.img_cache["3dpts"] = self.frames_dict[self.timestamp]["3dpts"]

            # Updates flag
            self.are_3d_points_loaded = True
            print("3d points loaded.")

    def get_keypts(self):
        """
        Calculates keypoints of RGB frame.
        """
        if "kps" not in list(self.frames_dict[self.timestamp].keys()):
            # Calculates keypoints
            kps = self.orb_detector.detect(self.img_cache["rgb"], None)

            # Filter keypoints
            kps = self.keypts_filter(kps)

            # Saves them onto dictionary of frames.
            self.frames_dict[self.timestamp]["kps"] = kps
            self.amount_img_calc += 1
            self.were_3dpts_calc = True
            print("Saved keypoints onto RAM.")

    def check_keypts(self):
        """
        Checks if key points are loaded onto image cache.
        """
        # Loads keypoints onto image cache
        if not self.are_kp_loaded:
            # Retrieves them from dictionary
            self.img_cache["kps"] = self.frames_dict[self.timestamp]["kps"]

            # Updates flag
            self.are_kp_loaded = True
            print("Keypoints loaded.")

    def keypts_filter(self, keypoints):
        """
        Filters keypoints based on distance from each other.
        """
        # Keypoints need to be sorted by strength in descending order
        # before feeding to SSC.
        keypoints = sorted(keypoints, key=lambda x: x.response, reverse=True)

        # Filter keypoints
        selected_keypoints = ssc(
            keypoints,
            self.main_yaml["filter_kps"]["ret_pts"],
            self.main_yaml["filter_kps"]["tol"],
            self.frame_size[0], self.frame_size[1]
        )
        return selected_keypoints

    def get_u_v_from_3dpts(self):
        """
        Estimates 2d points from 3d points.
        """
        if "3dpts_mtx" not in list(self.frames_dict[self.timestamp].keys()):
            print("Pasa aca?")
            print("Lista?", list(self.frames_dict[self.timestamp].keys()))
            # Iterates over each pixel of 3d points matrix
            for u in range(self.img_cache["3dpts"].shape[0]):
                for v in range(self.img_cache["3dpts"].shape[1]):
                    # Gets 3d points of (u, v) pair of depth camera
                    x_obj, y_obj, z_obj = self.img_cache["3dpts"][u, v]

                    # Estimates 2d point in RGB camera
                    u_obj = int(
                        np.around((x_obj * self.rgb_cam_mtx[0][0] / z_obj)
                                  + self.rgb_cam_mtx[0][2]))
                    v_obj = int(
                        np.around((y_obj * self.rgb_cam_mtx[1][1] / z_obj)
                                  + self.rgb_cam_mtx[1][2]))

                    # Filter points outside allowed values
                    if (480 > u_obj >= 0) and (640 > v_obj >= 0):
                        # Eliminates duplicated points
                        if (u_obj, v_obj) not in self.coord_values:
                            # Saves 3d points onto dictionary of coordinates
                            self.coord_values[
                                str(u_obj) + ", " + str(v_obj)] = (
                                    x_obj, y_obj, z_obj
                                )
                            debug_obj = False
                            if debug_obj:
                                if z_obj > 0:
                                    print("pos coord: ", u_obj, v_obj)
                                    print("real: ", x_obj, y_obj, z_obj)
            self.frames_dict[self.timestamp]["3dpts_mtx"] = self.coord_values
        else:
            print("Carga coordenadas")
            self.coord_values = self.frames_dict[self.timestamp]["3dpts_mtx"]

    def match_keypts_2dpts(self):
        """
        Matchs keypoints with 2d point estimations from 3d.
        """
        output_kp_list = []
        output_3dp_list = []
        # For each keypoint:
        for kp in self.img_cache["kps"]:
            # Gets (u, v) coordinates of keypoint
            u_obj, v_obj = int(kp.pt[0]), int(kp.pt[1])

            # If keypoints coordinates are in dictionary of 2d point
            # coordinates:
            if (str(u_obj) + ", " + str(v_obj)) in self.coord_values:
                # Retrieves 3D Coordinates on RGB frame
                str_x = self.coord_values[
                    str(u_obj) + ", " + str(v_obj)][0]
                str_y = self.coord_values[
                    str(u_obj) + ", " + str(v_obj)][1]

                # Fix Z axis scaling
                another_z = get_distance_meters(
                    self.img_cache["depth"][v_obj, u_obj])

                # Fix X and Y axis scaling
                if str_x < 0.2 and another_z < 2:
                    str_x = str_x * 2
                else:
                    str_x = str_x * 3
                str_y = str_y * 2.7
                str_x = np.float16(str_x)
                str_y = np.float16(str_y)
                another_z = np.float16(another_z)

                # Test another way of calculating Z:
                test_another_z = False
                if test_another_z:
                    another_z_2 = another_z * 480 / (2 * np.tan(45.6/2))
                    another_z_2 = another_z_2 / np.sqrt(
                        u_obj ** 2 + v_obj ** 2 + 480 / (
                            2 * np.tan(45.6/2))**2)
                    if another_z_2 < 0:
                        continue
                if another_z < 0:
                    continue

                test_some_vals = False
                if test_some_vals:
                    test_z_x = np.sqrt(another_z**2 - str_x**2)
                    test_z_y = np.sqrt(another_z**2 - str_y**2)
                    prom_z = np.mean((test_z_x, test_z_y))
                    print(prom_z)

                # Draw keypoints
                self.img_cache["rgb"] = draw_keypoints(
                    self.img_cache["rgb"],
                    [kp])

                # Show coordinates on screen
                show_text = False
                if show_text:
                    cv2.putText(
                        self.img_cache["rgb"],
                        # f"({str_x:.2f}, {str_y:.2f}, {str_z:.8f})",
                        f"({str_x:.2f}, {str_y:.2f}, {another_z:.2f})",
                        (int(u_obj) - 3, int(v_obj) - 3),
                        cv2.FONT_HERSHEY_SIMPLEX,
                        0.8, color=(0, 0, 255), thickness=2)
                output_kp_list.append(kp)
                output_3dp_list.append((str_x, str_y, another_z))
        self.frames_dict[self.timestamp]["kps"] = output_kp_list
        self.list_3dpts = output_3dp_list

    def invert_depth_image(self):
        """
        Inverts depth image.
        """
        self.img_cache["depth"] = cv2.bitwise_not(self.img_cache["depth"])

    def convert_depth_16b_2_8bit(self):
        """
        Converts depth image to seeable image.
        """
        self.img_cache["depth"] = self.img_cache["depth"].astype(np.uint8)

    def fill_misdetect_depth(self):
        """
        Fill missdetections depth.
        """
        depth_impaint = cv2.inpaint(self.img_cache["depth"],
                                    self.img_cache["mask"].astype(np.uint8), 2,
                                    cv2.INPAINT_TELEA)
        cv2.imshow("a", self.img_cache["mask"].astype(np.uint8))
        self.img_cache["depth"] = depth_impaint

    def threshold_mask(self):
        """
        Generates threshold mask
        """
        _, mask = cv2.threshold(
            self.img_cache["depth"], 2046, 2047, cv2.THRESH_BINARY)
        self.img_cache["mask"] = mask

    def generate_3D_array(self, key_press):
        """
        Generates array of 3D Points.
        """
        # Press "l" to begin generation process.
        if key_press == ord("l"):
            print("images processed: ", len(self.list_3dpts))
            if self.amount_img_calc != len(
                    self.main_yaml["frames"]["list_timestamps"]):
                print("Not enough processed images! Amount remaining: ",
                      len(self.main_yaml["frames"]["list_timestamps"])
                      - self.amount_img_calc)
            current_df = pd.DataFrame(self.list_3dpts, columns=["x", "y", "z"])
            current_df.to_csv(
                "data/timestamps/csv_files/"
                + "timestamp_" + self.timestamp + ".csv",
                index=False)

    def begin_measurement_process(self, key_press):
        """
        Begin measurement process.
        """
        if key_press == ord("v") and not self.flags["is recording"]:
            total_countdown = 5
            for idx in range(total_countdown):
                print('\007')
                print("Beginning recording in ", total_countdown - idx)
                time.sleep(1)
            print("Recording has started!")
            self.flags["is recording"] = True

    def measurement_process(self):
        """
        Run measurement process.
        """
        current_df = pd.DataFrame(self.list_3dpts, columns=["x", "y", "z"])
        current_df.to_csv(
            "data/rgbd/"
            + self.main_yaml["frames"]["chosen_id"] + "/csv_files/"
            + "timestamp_" + self.timestamp + ".csv",
            index=False)
        print(
            "data/rgbd/"
            + self.main_yaml["frames"]["chosen_id"] + "/csv_files/"
            + "timestamp_" + self.timestamp + ".csv saved!")
        # Updates current id, current timestamp, and image cache.
        (self.img_cache,
            self.current_id,
            self.timestamp) = change_img(
            1,
            self.img_cache,
            self.current_id,
            self.main_yaml["frames"]["list_timestamps"],
            self.frames_dict,
            self.timestamp)
        self.threshold_mask()
        # Resets flags so that 3d points and keypoints can be loaded.
        self.are_3d_points_loaded = False
        self.are_kp_loaded = False
        self.are_keypts_shown = False
        self.is_depth_processed = False
        self.were_3dpts_calc = False

    def show_frames(self):
        """
        Shows frames saved on image cache.
        """
        for frame_type in ["rgb", "depth"]:
            cv2.imshow(frame_type,
                       self.img_cache[frame_type])

    def take_capture(self, key_press):
        """
        Takes photo of current frame.
        """
        if key_press == ord("c"):
            current_time = datetime.now()
            timestamp_obj = current_time.strftime("%d_%b_%Y_%H_%M_%S")
            for str_type in ["rgb", "depth"]:
                # Saves frame with corresponding type.
                # E.g. frame_rgb_23847234.npy corresponds to a rgb image.
                cv2.imwrite(
                    self.main_yaml["frames"]["path_obj"] + "/"
                    + "frame_" + str_type + "_" + str(timestamp_obj) + ".png",
                    self.img_cache[str_type])
                print("Saved image ",
                      self.main_yaml["frames"]["path_obj"] + "/"
                      + "frame_" + str_type + "_" + str(timestamp_obj))

    def close_program(self, key_press):
        """
        Closes program if key is pressed.
        """
        # Press "q" to quit process
        if key_press == ord("q"):
            cv2.destroyAllWindows()
            sys.exit(1)

    def convert_depth(self):
        """
        Converts raw depth values to meters.
        """
        for u in range(self.img_cache["depth"].shape[0]):
            for v in range(self.img_cache["depth"].shape[1]):
                self.img_cache["depth"][u, v] = get_distance_meters(
                    self.img_cache["depth"][u, v])

    def save_3dpoints_mtx(self):
        """
        Given key pressed, asks for input data so it can interpolate movement.
        """
        # if key_press == ord("j"):
        current_dict = self.frames_dict[self.timestamp]["3dpts_mtx"]
        json_obj = json.dumps(current_dict)
        file_obj = open(
            "data/rgbd/"
            + self.main_yaml["frames"]["chosen_id"] + "/3dpts_mtx/"
            + "timestamp_" + self.timestamp + ".json",
            "w", encoding="utf-8"
        )
        file_obj.write(json_obj)
        print(
            "data/rgbd/"
            + self.main_yaml["frames"]["chosen_id"] + "/3dpts_mtx/"
            + "timestamp_" + self.timestamp + ".json saved!")
        file_obj.close()
        # Updates current id, current timestamp, and image cache.
        (self.img_cache,
            self.current_id,
            self.timestamp) = change_img(
            1,
            self.img_cache,
            self.current_id,
            self.main_yaml["frames"]["list_timestamps"],
            self.frames_dict,
            self.timestamp)
        self.threshold_mask()
        # Resets flags so that 3d points and keypoints can be loaded.
        self.are_3d_points_loaded = False
        self.are_kp_loaded = False
        self.are_keypts_shown = False
        self.is_depth_processed = False
        self.were_3dpts_calc = False

    def input_pose_estimation(self, key_press):
        """
        Inputs pose estimation for current frame.
        """
        if key_press == ord("u"):
            list_pose_obj = []
            while len(list_pose_obj) != 7:
                list_pose_obj = input(
                    "Ingrese valores para x, y, z, i, j, k, w")
                list_pose_obj = list_pose_obj.split()
                list_pose_obj = [float(item_obj) for item_obj in list_pose_obj]
                if len(list_pose_obj) != 7:
                    print("No se tienen 7 valores. Intente nuevamente")
            self.df_pose.loc[self.current_id] = list_pose_obj
            self.df_pose.to_csv(
                "data/rgbd/" + self.main_yaml["frames"]["chosen_id"]
                + "/df_pose_" + self.main_yaml["frames"]["chosen_id"] + ".csv",
                index=False)

    def load_3dpts_mtx(self):
        """
        Checks if 3dpts_mtx file exists. If it is True, loads it and stores it
        onto self.coord_values.
        """
        for idx, timestamp_obj in enumerate(
                self.main_yaml["frames"]["list_timestamps"]):
            if idx == 300:
                break
            folder_obj = ("/home/diego/RFS_SLAM/rfsslam/"
                          + "data/rgbd/"
                          + self.main_yaml["frames"]["chosen_id"]
                          + "/3dpts_mtx/")
            bool_file_exists = os.path.exists(
                folder_obj + "timestamp_" + timestamp_obj + ".json")
            if bool_file_exists:
                file_obj = open(
                    folder_obj + "timestamp_" + timestamp_obj + ".json",
                    encoding="utf-8")
                dict_obj = json.load(file_obj)
                self.frames_dict[timestamp_obj]["3dpts_mtx"] = dict_obj
                file_obj.close()
                print("Loaded file n", idx)
            else:
                print("File n ", idx, " doesn't exist!")

    def get_list_frames_from_folder(self):
        """
        Returns list of all frames on "chosen_id" in yaml.
        """
        main_path = "data/rgbd/"
        folder_ids = os.listdir(main_path)
        if self.main_yaml["frames"]["chosen_id"] in folder_ids:
            frames_list = os.listdir(
                main_path + self.main_yaml["frames"]["chosen_id"]
                + "/images/rgb/")
            frames_list = [name_obj[10:-4] for name_obj in frames_list]
            print("N frames loaded:", len(frames_list))
            frames_list.sort()
            print("sorted list: ", frames_list)
            self.main_yaml["frames"]["list_timestamps"] = frames_list
        else:
            sys.exit(
                "Chosen id: " + str(self.main_yaml["frames"]["chosen_id"])
                + "not in data/rgbd/.")

    def copy_current_file_to_new_timestamp(self, key_press):
        """
        Copy current timestamp to folder: new_timestamp.
        """
        if key_press == ord("t"):
            src_folder = (
                "data/rgbd/" + self.main_yaml["frames"]["chosen_id"]
                + "/images/")
            dst_folder = (
                "data/rgbd/" + self.new_timestamp + "/images/")

            for frame_type in ["rgb", "depth", "ir"]:
                # file names
                src_file = (
                    src_folder + frame_type + "/" + "frame_" + frame_type
                    + "_" + self.timestamp + ".npy")
                dst_file = (
                    dst_folder + frame_type + "/" + "frame_" + frame_type
                    + "_" + self.timestamp + ".npy")

                try:
                    # copy file
                    shutil.copyfile(src_file, dst_file)
                    # destination folder after copying
                    print("Destination after copying", os.listdir(dst_folder))
                except SameFileError:
                    print("We are trying to copy the same File")
                except IsADirectoryError:
                    print("The destination is a directory")