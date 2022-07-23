"""
main_script.py:
    Main script.
"""
from utils.misc_tools import (
    open_calibr_dict, open_cfg_yaml
)
from classes.frames_class import FrameDisplay
from classes.realt_class import VideoDisplay
from classes.class_calibr import CalibrDisplay


def main():
    """
    Executes main code.
    """
    # Initialize YAML configuration dictionary
    main_yaml = open_cfg_yaml()

    # Get dict containing RGB, depth and stereo parameters
    dict_calibr = open_calibr_dict()

    # Main loop
    while True:

        # Frame by frame mode.
        if main_yaml["mode"] == "framebyframe":

            frame_obj = FrameDisplay(main_yaml, dict_calibr)
            frame_obj.load_first_frame()
            frame_obj.run_algorithms()

        # Real time mode
        if main_yaml["mode"] == "realtime":
            frame_obj = VideoDisplay(main_yaml)
            frame_obj.run_algorithms()

        # Calibration mode.
        if main_yaml["mode"] == "calibration":

            frame_obj = CalibrDisplay(main_yaml)
            frame_obj.run_algorithms()


if __name__ == "__main__":
    main()
