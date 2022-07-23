"""
A
"""
import json
import numpy as np
import yaml


def open_calibr_dict():
    """
    Opens calibration dictionary.
    """
    # Path with camera YAML.
    path_obj = "src/python_script/config/camera_calibration_params.JSON"
    with open(path_obj, encoding='utf8') as json_obj:
        output_dict = json.load(json_obj)
    return output_dict


def open_cfg_yaml():
    """
    Opens configuration YAML.

    Parameters
    ----------
    None.

    Returns
    -------
    output_yaml:
        Output YAML.
    """
    yaml_name = "cfg_parameters.yaml"
    with open(
            "src/python_script/config/" + yaml_name, "r", encoding='utf8'
            ) as yml:
        output_yaml = yaml.load(
            yml, Loader=yaml.FullLoader)
    return output_yaml


def open_yaml_cameras(str_obj):
    """
    Opens YAML file of RGB/depth camera.

    Parameters
    ----------
    str_obj: str
        Camera to retrieve parameters. Can be "rgb" or "depth".

    Returns
    -------
    output_dict: dict
        Dictionary with camera parameters.
    """
    # Path with camera YAML.
    path_obj = ('/home/diego/.ros/camera_info/'
                + str_obj + '_B00365203483122B.yaml')
    with open(path_obj, encoding='utf8') as file:
        output_dict = yaml.load(file, Loader=yaml.FullLoader)
    return output_dict


def save_frame(path_obj, frame_obj, str_type, timestamp_obj, main_yaml):
    """
    Saves frame of image onto disk.

    Parameters
    ----------
    frame_obj: np.ndarray
        Frame to save.
    str_type: str
        Type of frame to save.
    timestamp_obj: float
        Timestamp.
    main_yaml:
    """
    # Saves frame with corresponding type.
    # E.g. frame_rgb_23847234.npy corresponds to a rgb image.
    np.save(path_obj + "/"
            + "frame_" + str_type + "_" + str(timestamp_obj) + ".npy",
            frame_obj)

    if str(timestamp_obj) not in main_yaml["frames"]["list_timestamps"]:
        main_yaml["frames"]["list_timestamps"].append(str(timestamp_obj))
        with open(
                "src/python_script/config/cfg_parameters.yaml",
                "w", encoding="utf-8"
                ) as yaml_file:
            yaml_file.write(yaml.dump(main_yaml, default_flow_style=False))
    return main_yaml


def load_frames(path_obj, list_timestamps):
    """
    Loads frame given path and name of the file.

    Parameters
    ----------
    path_obj: str
        Main path of folder that contains frames.
    list_timestamps: list
        List with available timestamps.

    Returns
    -------
    dict_frames: dict
        Dictionary with frames saved given timestamp.
        Format is: dict[timestamps][type_frame]
    """
    dict_frames = {}
    # Iterates over all saved timestamps
    for timestamp_obj in list_timestamps:
        dict_frames[str(timestamp_obj)] = {}
        # Saves all types of frames given timestamp
        for str_type in ["rgb", "depth", "ir"]:
            dict_frames[str(timestamp_obj)][str_type] = np.load(
                path_obj + str_type + "/" + "frame_" + str_type + "_"
                + str(timestamp_obj) + ".npy"
            )
    return dict_frames


def change_img(int_direction, img_cache, current_id, list_folder, frames_dict,
               current_timestamp):
    """
    Rotates over timestamps in image cache.
    """
    # Modifies current ID with input.
    current_id += int_direction

    # If ID is bigger than ID of last image in the list:
    # Loops to the beggining
    if current_id == len(list_folder):
        current_id = 0
    elif current_id == -1:
        current_id = len(list_folder) - 1
    print("Current image index:", current_id + 1, "/", len(list_folder))
    current_timestamp = list_folder[current_id]
    # Loads image given ID
    for frame_type in ["rgb", "depth", "ir"]:
        img_cache[frame_type] = (
            frames_dict[list_folder[current_id]][frame_type])

    # Only to use to compare between two images.
    if current_id + 1 == len(list_folder):
        img_cache["nextrgb"] = (
            frames_dict[list_folder[0]]["rgb"])
    else:
        img_cache["nextrgb"] = (
            frames_dict[list_folder[current_id + 1]]["rgb"])
    return img_cache, current_id, current_timestamp


def merge_R_t_one_mtx(rot_mtx: np.ndarray,
                      trans_mtx: np.ndarray) -> np.ndarray:
    """
    Merges rotation vector and translation matrix.

    Parameters
    ----------
    rot_mtx : np.ndarray
        Rotation matrix.
    trans_mtx: np.ndarray
        Translation matrix.

    Returns
    -------
    output_mtx: np.ndarray
        Matrix with [R|t].
    """
    # Initializes output matrix
    output_mtx = np.zeros((3, 4))
    output_mtx[0:3, 0:3] = rot_mtx
    output_mtx[:, -1] = trans_mtx
    return output_mtx


class MyEncoder(json.JSONEncoder):
    """
    Encoder to convert ndarrays to list.
    """
    def default(self, o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        return json.JSONEncoder.default(self, o)
