Python scripts for database creation and processing
============================

## Introduction

This subfolder contains a series of scripts that handles the following tasks:
1. RGB-D Camera calibration.
2. Captures image/video input from camera.
3. ORB feature extraction and 3D estimation given input.
4. Measurement and pose dump onto csv.

Main file is called `main_script.py`. To configure it, go to `config/cfg_parameters.yaml`. There you can adjust mode (Camera calibration, camera input recording, or frame processing), assign database to analyze, and more (detailed description will be in its main section).

## Folders structure

## Requirements

Under construction.

## RGB-D Camera calibration

Under construction.

## Camera input recording

If you want to record a database from the RGB-D input, you need to set `mode` in `config/cfg_parameters.yaml` to **realtime**:
```
(...)
key_press:
  '110': 1
  '98': -1
mode: realtime
new_timestamp: 24_jul_2022_16_52_00
(...)
```

After that, launch the main script with the following command:

```
python3 main_script.py
```

If you have the RGB-D camera well configurated, it will launch the system right away, showing to you the following windows:
* **rgb**: Shows input from RGB camera
* **depth**: Shows input from depth camera
* **ir**: Shows input from infrared camera

You can move freely the camera, but, because of lack of optimization, good framerate is not assured (Kinect V1 registeres 1~2 fps).

You have **two modes of information storage**: by capturing photos pressing a button, and capturing a stream of photos.

This mode has a series of key presses that allows you a set of actions:
* **c**: Captures a *photo* of current input. By default, a folder with the date and time of the launch of the script is created, and these photos will be stored there.
    * The information will be separated in folders, given the nature of the inputs.
    * E.g: Suppose you run the script at 4/08/22, at 12:00:00. It will create a folder with the following name: `04_Aug_22_12_00_00`, which will store all the photos that you will take with the **c** key press.
    * Each time you press **c**, it will store in `04_Aug_22_12_00_00/rgb/` all rgb photos, `04_Aug_22_12_00_00/depth/` all depth photos, and `04_Aug_22_12_00_00/ir/` all infrared photos.
    * All photos will be saved with **.npy** format.

* **v**: Starts **input stream recording** in 5 seconds. It will save input stream of photos as fast as it can, storing them in a folder **with the date and time of the BEGINNING OF THE RECORDING**, **NOT THE DATE AND TIME OF THE SCRIPT LAUNCH!**
    * If you press **v** at 4/08/22, 12:10:10, it'll create a folder with the name of `04_Aug_22_12_10_15` that will contain all the information of this recording.
* **s**: **Stops recording mode**.
* **h**: Displays help commands (at 04/08 is not well implemented yet :c)
* **q**: Quits program.

## Frame processing

This mode **REQUIRES** you to have a database recorded, shall it be with screen capture or stream input recording. For you to activate this mode, you need to set both `mode` in `config/cfg_parameters.yaml` to **framebyframe** and also `chosen_id` to your folder database name.

E.g, if you will use folder `24_jul_2022_16_59_00`:
```
(...)
frames:
  chosen_id: 24_jul_2022_16_59_00
(...)
key_press:
  '110': 1
  '98': -1
mode: framebyframe
(...)
```

This mode mainly handles 3d coordinate extraction from keypoints calculated in certain folder.