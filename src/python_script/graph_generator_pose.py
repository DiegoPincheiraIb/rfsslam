"""
A
"""
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def delete_unwanted_cols(df_obj, list_cols):
    """
    Deletes unwanted columns.
    """
    for col_name in list_cols:
        del df_obj[col_name]


def pose_df_proccesing(df_obj):
    """
    Process pose dataframe.

    Parameters
    ----------
    df_obj: pd.DataFrame
        Pose dataframe.

    Returns
    -------
    filtered_df_obj: pd.DataFrame
        Processed pose dataframe.
    """
    df_obj["timestamp"] = df_obj["timestamp"].div(3)
    df_obj["weight"].fillna(1, inplace=True)
    df_obj = df_obj.drop_duplicates()
    lista_tstmp = list(collections.Counter(df_obj["timestamp"]))
    pose_filtered_data = []
    for timestamp in lista_tstmp:
        timestamp = int(timestamp)
        tstmp3 = (
            df_obj.loc[df_obj.loc[:, "timestamp"] == lista_tstmp[timestamp]])
        lolazo = tstmp3[tstmp3.weight == tstmp3.weight.max()]
        media_obj = list(lolazo.mean())[1:-1]
        if np.isnan(np.array(media_obj)).any():
            print("timestamp: ", timestamp)
            print("aqui: ", list(df_obj.iloc[timestamp, :]))
        pose_filtered_data.append(media_obj)
    filtered_df_obj = pd.DataFrame(
        pose_filtered_data, columns=["x", "y", "z", "qy", "qw"])
    return filtered_df_obj


# Read files
df_pose = pd.read_table(
    # "/home/diego/RFS_SLAM/rfsslam/build/data/rbphdslam6d/particlePose.dat",
    "/home/diego/RFS_SLAM/rbphdslam6d_true/particlePose.dat",
    sep="   ",
    names=[
        "timestamp", "particle",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])

df_GT = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/gtPose.dat",
    names=["timestamp", "x", "y", "z", "qx", "qy", "qz", "qw"], sep="   ")

df_odometry = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/deadReckoning.dat",
    names=["timestamp", "x", "y", "z", "qx", "qy", "qz", "qw"], sep="   ")

# -------------------------
folder_obj_1 = "rbphdslam6d_2Y"
cov_value = "2"
df_pose_1 = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/"
    + folder_obj_1 + "/particlePose.dat",
    sep="   ",
    names=[
        "timestamp", "particle",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])
df_odometry_1 = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/"
    + folder_obj_1 + "/deadReckoning.dat",
    sep="   ",
    names=[
        "timestamp",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])

folder_obj_2 = "rbphdslam6d_20X"
df_pose_2 = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/"
    + folder_obj_2 + "/particlePose.dat",
    sep="   ",
    names=[
        "timestamp", "particle",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])
df_odometry_2 = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/"
    + folder_obj_2 + "/deadReckoning.dat",
    sep="   ",
    names=[
        "timestamp",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])
print(df_odometry_1.head())
# ---------------UNUSED -------
df_pose_3 = pd.read_table(
    "/home/diego/RFS_SLAM/rbphdslam6d_true/rbphdslam6d_1200/particlePose.dat",
    sep="   ",
    names=[
        "timestamp", "particle",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])

delete_unwanted_cols(df_pose, ["particle", "qx", "qz"])
delete_unwanted_cols(df_pose_1, ["particle", "qx", "qz"])
delete_unwanted_cols(df_pose_2, ["particle", "qx", "qz"])
delete_unwanted_cols(df_pose_3, ["particle", "qx", "qz"])
delete_unwanted_cols(df_GT, ["timestamp", "qx", "qz"])
delete_unwanted_cols(df_odometry, ["timestamp", "qx", "qz"])
delete_unwanted_cols(df_odometry_1, ["timestamp", "qx", "qz"])
delete_unwanted_cols(df_odometry_2, ["timestamp", "qx", "qz"])

filtered_df_pose = pose_df_proccesing(df_pose)
filtered_df_pose_1 = pose_df_proccesing(df_pose_1)
filtered_df_pose_2 = pose_df_proccesing(df_pose_2)
filtered_df_pose_3 = pose_df_proccesing(df_pose_3)

filtered_df_pose.qw[0] = 1
filtered_df_pose_1.qw[0] = 1
filtered_df_pose_2.qw[0] = 1

error_gt_odo = df_GT.subtract(df_odometry).abs()
error_gt_odo_1 = df_GT.subtract(df_odometry_1).abs()
error_gt_odo_2 = df_GT.subtract(df_odometry_2).abs()
error_gt_pose = df_GT.subtract(filtered_df_pose).abs()
error_gt_pose_1 = df_GT.subtract(filtered_df_pose_1).abs()
error_gt_pose_2 = df_GT.subtract(filtered_df_pose_2).abs()
error_gt_pose_3 = df_GT.subtract(filtered_df_pose_3).abs()
lista_idxs = list(error_gt_pose.index)
asd = False
if asd:
    euclidean_odo = np.sqrt(
        np.square(df_GT.subtract(df_odometry).x)
        + np.square(df_GT.subtract(df_odometry).y)
        + np.square(df_GT.subtract(df_odometry).z))

    euclidean_est = np.sqrt(
        np.square(df_GT.subtract(filtered_df_pose).x)
        + np.square(df_GT.subtract(filtered_df_pose).y)
        + np.square(df_GT.subtract(filtered_df_pose).z))

# Print errores
print_info_filtered_df_pose = False
if print_info_filtered_df_pose:
    print(filtered_df_pose)
    print(df_odometry)
    print(error_gt_pose["x"].isnull().sum())

# Convertir valores de odometría a listas para poder graficarlos
x_values_odometry = list(df_odometry["x"])
x_values_odometry = [value * -1 for value in x_values_odometry]
x_values_odometry_1 = list(df_odometry_1["x"])
x_values_odometry_1 = [value * -1 for value in x_values_odometry_1]
x_values_odometry_2 = list(df_odometry_2["x"])
x_values_odometry_2 = [value * -1 for value in x_values_odometry_2]
y_values_odometry = list(df_odometry["y"])
z_values_odometry = list(df_odometry["z"])
z_values_odometry_1 = list(df_odometry_1["z"])
z_values_odometry_2 = list(df_odometry_2["z"])

x_values_GT = list(df_GT["x"])
y_values_GT = list(df_GT["y"])
z_values_GT = list(df_GT["z"])

x_values_pose = list(filtered_df_pose["x"])
x_values_pose_1 = list(filtered_df_pose_1["x"])
x_values_pose_2 = list(filtered_df_pose_2["x"])

y_values_pose = list(filtered_df_pose["y"])
y_values_pose_1 = list(filtered_df_pose_1["y"])
y_values_pose_2 = list(filtered_df_pose_2["y"])

z_values_pose = list(filtered_df_pose["z"])
z_values_pose_1 = list(filtered_df_pose_1["z"])
z_values_pose_2 = list(filtered_df_pose_2["z"])

x_values_GT = [value * -1 for value in x_values_GT]
x_values_pose = [value * -1 for value in x_values_pose]
x_values_pose_1 = [value * -1 for value in x_values_pose_1]
x_values_pose_2 = [value * -1 for value in x_values_pose_2]
# Frame donde se retorna a la posición inicial
frame_beg = 176

print_console_error = True
if print_console_error:
    print(
        "error en ejes: ",
        error_gt_pose.loc[frame_beg, "x"],
        error_gt_pose.loc[frame_beg, "y"],
        error_gt_pose.loc[frame_beg, "z"],
        error_gt_pose.loc[frame_beg, "qy"],
        error_gt_pose.loc[frame_beg, "qw"])
    print(
        "Coord en ejes de odometria: ",
        df_odometry.loc[frame_beg, "x"],
        df_odometry.loc[frame_beg, "y"],
        df_odometry.loc[frame_beg, "z"],
        df_odometry.loc[frame_beg, "qy"],
        df_odometry.loc[frame_beg, "qw"])
    print(
        "Coord estimadas: ",
        filtered_df_pose.loc[frame_beg, "x"],
        filtered_df_pose.loc[frame_beg, "y"],
        filtered_df_pose.loc[frame_beg, "z"],
        filtered_df_pose.loc[frame_beg, "qy"],
        filtered_df_pose.loc[frame_beg, "qw"])


plot_dict = {
    "Error comparacion quat": False,
    "Error comparacion cov": False,
    "Error raw cart": False,
    "Error raw quat": False,
    "Error acumulado x": False,
    "Error acumulado y": False,
    "Error acumulado z": False,
    "Trayectoria": False,
    "TrayectoriaX": False,
    "TrayectoriaY": True,
    "TrayectoriaZ": False,
}

# Plot graphics
plot_graphs = True
if plot_graphs:
    if plot_dict["Error comparacion quat"]:
        axis_obj = "Qw"
        axis_obj_minus = "qw"
        fig_obj, ax_obj = plt.subplots(2, 1, sharex='col', sharey='row')
        ax_obj[0].plot(
            lista_idxs,
            error_gt_odo.loc[:, axis_obj_minus].cumsum(), "--",
            label=("Error acumulado en " + axis_obj + " para odometría"),
            color="green")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_odo_1.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en " + axis_obj
                + " para 1.5" + "\u03C3 en odometría."),
            color="darkgreen")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_odo_2.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en " + axis_obj
                + " para 2" + "\u03C3 en odometría."),
            color="springgreen")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_pose.loc[:, axis_obj_minus].cumsum(), "--",
            label=(
                "Error acumulado en " + axis_obj
                + " para " + "\u03C3" + " original"),
            color="red")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_pose_1.loc[:, axis_obj_minus].cumsum(),
            label="Error acumulado en " + axis_obj + " para 1.5" + "\u03C3",
            color="darkblue")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_pose_2.loc[:, axis_obj_minus].cumsum(),
            label="Error acumulado en " + axis_obj + " para 2" + "\u03C3",
            color="orange")
        # Plot points of return
        ax_obj[0].plot(
            frame_beg,
            error_gt_pose.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko",
            label="Retorno a la posición inicial")
        ax_obj[0].plot(
            frame_beg,
            error_gt_pose_1.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[0].plot(
            frame_beg,
            error_gt_pose_2.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[0].plot(
            frame_beg,
            error_gt_odo.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[0].plot(
            frame_beg,
            error_gt_odo_1.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[0].plot(
            frame_beg,
            error_gt_odo_2.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[0].set_ylabel("Error acumulado en Qw[m]")
        ax_obj[0].legend(loc='best')
        ax_obj[0].set_title(
            "Error acumulado en los cuaterniones Qw y Qy a lo largo de"
            + "\n la trayectoria para diferentes valores de covarianzas.")
        # plt.subplot(312)
        axis_obj = "Qy"
        axis_obj_minus = "qy"
        ax_obj[1].plot(
            lista_idxs,
            error_gt_odo.loc[:, axis_obj_minus].cumsum(), "--",
            label=("Error acumulado en " + axis_obj + " para odometría"),
            color="green")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_odo_1.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en " + axis_obj
                + " para 1.5" + "\u03C3 en odometría."),
            color="darkgreen")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_odo_2.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en " + axis_obj
                + " para 2" + "\u03C3 en odometría."),
            color="springgreen")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_pose.loc[:, axis_obj_minus].cumsum(), "--",
            label=(
                "Error acumulado en " + axis_obj
                + " para " + "\u03C3" + " original"),
            color="red")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_pose_1.loc[:, axis_obj_minus].cumsum(),
            label="Error acumulado en " + axis_obj + " para 1.5" + "\u03C3",
            color="darkblue")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_pose_2.loc[:, axis_obj_minus].cumsum(),
            label="Error acumulado en " + axis_obj + " para 2" + "\u03C3",
            color="orange")
        # Plot points of return
        ax_obj[1].plot(
            frame_beg,
            error_gt_pose.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko",
            label="Retorno a la posición inicial")
        ax_obj[1].plot(
            frame_beg,
            error_gt_pose_1.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[1].plot(
            frame_beg,
            error_gt_pose_2.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[1].plot(
            frame_beg,
            error_gt_odo.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[1].plot(
            frame_beg,
            error_gt_odo_1.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        ax_obj[1].plot(
            frame_beg,
            error_gt_odo_2.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        plt.subplots_adjust(hspace=.0)
        ax_obj[1].set_ylabel("Error en Qy[m]")
        ax_obj[1].legend(loc='best')
        plt.xlabel("Tiempo [frames]")
        plt.ylabel("Error acumulado en eje " + axis_obj + " [m]")

    if plot_dict["Error comparacion cov"]:
        axis_obj = "Y"
        axis_obj_minus = "y"

        # print("Son iguales?", error_gt_odo_1.x.equals(error_gt_odo_2.x))
        plt.figure(figsize=(7, 4))
        plt.plot(
            lista_idxs, error_gt_odo.loc[:, axis_obj_minus].cumsum(), "--",
            label="Error acumulado en eje " + axis_obj + " para odometría",
            color="green")
        plt.plot(
            lista_idxs, error_gt_odo_1.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en eje " + axis_obj
                + " para 1.5" + "\u03C3 en odometría."),
            color="darkgreen")
        plt.plot(
            lista_idxs, error_gt_odo_2.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en eje " + axis_obj
                + " para 2" + "\u03C3 en odometría."),
            color="springgreen")
        plt.plot(
            lista_idxs, error_gt_pose.loc[:, axis_obj_minus].cumsum(), "--",
            label=(
                "Error acumulado en eje " + axis_obj
                + " para " + "\u03C3" + " original"),
            color="red")
        plt.plot(
            lista_idxs, error_gt_pose_1.loc[:, axis_obj_minus].cumsum(),
            label=(
                "Error acumulado en eje "
                + axis_obj + " para 1.5" + "\u03C3"),
            color="darkblue")
        plt.plot(
            lista_idxs, error_gt_pose_2.loc[:, axis_obj_minus].cumsum(),
            label="Error acumulado en eje " + axis_obj + " para 2" + "\u03C3",
            color="orange")

        # Plot points of return
        plt.plot(
            frame_beg,
            error_gt_pose.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko",
            label="Retorno a la posición inicial")
        plt.plot(
            frame_beg,
            error_gt_pose_1.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        plt.plot(
            frame_beg,
            error_gt_pose_2.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        plt.plot(
            frame_beg,
            error_gt_odo.loc[:, axis_obj_minus].cumsum()[frame_beg], "ko")
        plt.legend(loc='upper left')
        plt.xlabel("Tiempo [frames]")
        plt.ylabel("Error acumulado en eje " + axis_obj + " [m]")
        plt.title(
            "Error acumulado en el eje " + axis_obj + " a lo largo de"
            + "\n la trayectoria para diferentes valores de covarianza.")

    if plot_dict["Error raw cart"]:
        fig_obj, ax_obj = plt.subplots(3, 1, sharex='col', sharey='row')
        ax_obj[0].plot(
            lista_idxs, [0]*len(lista_idxs), label="Referencia")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_odo.x, label="Error de dead reckoning", color="green")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_pose.x,
            label="Error de trayectoria estimada", color="red")
        ax_obj[0].plot(
            frame_beg, 0, "ko",
            label="Retorno a la posición inicial")
        ax_obj[0].set_ylabel("Error en eje X[m]")
        ax_obj[0].legend(loc='best')
        ax_obj[0].set_title("Error estimado a lo largo de la trayectoria")
        # plt.subplot(312)
        ax_obj[1].plot(
            lista_idxs, [0]*len(lista_idxs), label="Referencia")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_odo.y, label="Error de dead reckoning", color="green")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_pose.y,
            label="Error de trayectoria estimada", color="red")
        ax_obj[1].set_ylabel("Error en eje Y[m]")
        ax_obj[1].plot(
            frame_beg, 0, "ko",
            label="Retorno a la posición inicial")
        ax_obj[1].legend(loc='best')
        # plt.subplot(313)
        ax_obj[2].plot(
            lista_idxs, [0]*len(lista_idxs), label="Referencia")
        ax_obj[2].plot(
            lista_idxs, error_gt_odo.z,
            label="Error de dead reckoning", color="green")
        ax_obj[2].plot(
            lista_idxs, error_gt_pose.z,
            label="Error de trayectoria estimada", color="red")
        plt.subplots_adjust(hspace=.0)
        ax_obj[2].plot(
            frame_beg, 0, "ko",
            label="Retorno a la posición inicial")
        ax_obj[2].legend(loc='best')
        plt.ylabel("Error en eje Z[m]")
        plt.xlabel("Tiempo [frames]")

    if plot_dict["Error raw quat"]:
        fig_obj, ax_obj = plt.subplots(2, 1, sharex='col', sharey='row')
        ax_obj[0].plot(
            lista_idxs, [0]*len(lista_idxs), label="Referencia")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_odo.qw, label="Error de dead reckoning", color="green")
        ax_obj[0].plot(
            lista_idxs,
            error_gt_pose.qw,
            label="Error de trayectoria estimada", color="red")
        ax_obj[0].plot(
            frame_beg, 0, "ko", label="Retorno a la posición inicial")
        ax_obj[0].set_ylabel("Error en Qw[m]")
        ax_obj[0].legend(loc='best')
        ax_obj[0].set_title(
            "Error estimado a lo largo de la trayectoria en cuaterniones")
        # plt.subplot(312)
        ax_obj[1].plot(
            lista_idxs, [0]*len(lista_idxs), label="Referencia")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_odo.qy, label="Error de dead reckoning", color="green")
        ax_obj[1].plot(
            lista_idxs,
            error_gt_pose.qy,
            label="Error de trayectoria estimada", color="red")
        ax_obj[1].set_ylabel("Error en Qy[m]")
        ax_obj[1].plot(
            frame_beg, 0, "ko", label="Retorno a la posición inicial")
        ax_obj[1].legend(loc='best')
        plt.ylabel("Error en Qy[m]")
        plt.xlabel("Tiempo [frames]")

    if plot_dict["Error acumulado x"]:
        plt.figure(figsize=(7, 4))
        plt.plot(
            lista_idxs, error_gt_odo.x.cumsum(), "--",
            label="Error acumulado en eje X para odometría", color="green")
        plt.plot(
            lista_idxs, error_gt_pose_2.x.cumsum(),
            label="Error acumulado en eje X para p=20", color="darkblue")
        plt.plot(
            lista_idxs, error_gt_pose_1.x.cumsum(),
            label="Error acumulado en eje X para p=200", color="orange")
        plt.plot(
            lista_idxs, error_gt_pose.x.cumsum(), "--",
            label="Error acumulado en eje X para p=600", color="red")
        plt.plot(
            frame_beg, error_gt_pose.x.cumsum()[frame_beg], "ko",
            label="Retorno a la posición inicial")
        plt.plot(frame_beg, error_gt_pose_1.x.cumsum()[frame_beg], "ko")

        plt.plot(frame_beg, error_gt_pose_2.x.cumsum()[frame_beg], "ko")
        plt.plot(frame_beg, error_gt_odo.x.cumsum()[frame_beg], "ko")
        plt.legend(loc='upper left')
        plt.xlabel("Tiempo [frames]")
        plt.ylabel("Error acumulado en eje X [m]")
        plt.title(
            "Error acumulado en el eje X a lo largo de"
            + "\n la trayectoria para diferente cantidad de partículas.")

    if plot_dict["Error acumulado y"]:
        plt.figure(figsize=(7, 4))
        plt.plot(
            lista_idxs, error_gt_odo.y.cumsum(), "--",
            label="Error acumulado en eje Y para odometría", color="green")
        plt.plot(
            lista_idxs, error_gt_pose_1.y.cumsum(),
            label="Error acumulado en eje Y para p=20", color="darkblue")
        plt.plot(
            lista_idxs, error_gt_pose_2.y.cumsum(),
            label="Error acumulado en eje Y para p=200", color="orange")
        plt.plot(
            lista_idxs, error_gt_pose.y.cumsum(), "--",
            label="Error acumulado en eje Y para p=600", color="red")
        plt.plot(
            frame_beg, error_gt_pose.y.cumsum()[frame_beg], "ko",
            label="Retorno a la posición inicial")
        plt.plot(frame_beg, error_gt_pose_1.y.cumsum()[frame_beg], "ko")

        plt.plot(frame_beg, error_gt_pose_2.y.cumsum()[frame_beg], "ko")
        plt.plot(frame_beg, error_gt_odo.y.cumsum()[frame_beg], "ko")
        plt.legend(loc='upper left')
        plt.xlabel("Tiempo [frames]")
        plt.ylabel("Error acumulado en eje Y [m]")
        plt.title(
            "Error acumulado en el eje Y a lo largo de"
            + "\n la trayectoria para diferente cantidad de partículas.")

    if plot_dict["Error acumulado z"]:
        plt.figure(figsize=(7, 4))
        plt.plot(
            lista_idxs, error_gt_odo.z.cumsum(), "--",
            label="Error acumulado en eje Z para odometría", color="green")
        plt.plot(
            lista_idxs, error_gt_pose_2.z.cumsum(),
            label="Error acumulado en eje Z para p=20", color="darkblue")
        plt.plot(
            lista_idxs, error_gt_pose_1.z.cumsum(),
            label="Error acumulado en eje Z para p=200", color="orange")
        plt.plot(
            lista_idxs, error_gt_pose.z.cumsum(), "--",
            label="Error acumulado en eje Z para p=600", color="red")
        plt.plot(
            frame_beg, error_gt_pose.z.cumsum()[frame_beg], "ko",
            label="Retorno a la posición inicial")
        plt.plot(frame_beg, error_gt_pose_1.z.cumsum()[frame_beg], "ko")
        plt.plot(frame_beg, error_gt_pose_2.z.cumsum()[frame_beg], "ko")
        plt.plot(frame_beg, error_gt_odo.z.cumsum()[frame_beg], "ko")
        plt.legend(loc='upper left')
        plt.xlabel("Tiempo [frames]")
        plt.ylabel("Error acumulado en eje Z [m]")
        plt.title(
            "Error acumulado en el eje Z a lo largo de \n"
            + "la trayectoria para diferente cantidad de partículas.")

    if plot_dict["Trayectoria"]:
        plt.figure(figsize=(6.5, 7))
        plt.plot(
            x_values_GT, z_values_GT, label="Ground truth")
        plt.plot(
            x_values_odometry, z_values_odometry, "--",
            label=("Dead reckoning para " + "\u03C3" + " original en Z"),
            color="green")
        plt.plot(
            x_values_odometry_1, z_values_odometry_1,
            label=(
                "Dead reckoning" + " para " + cov_value + "\u03C3" + " en Z"),
            color="darkgreen")
        plt.plot(
            x_values_pose, z_values_pose, "--",
            label=(
                "Pose estimada del vehículo para " + "\u03C3"
                + " original en Z"),
            color="red")
        plt.plot(
            x_values_pose_1, z_values_pose_1,
            label=(
                "Pose estimada del vehículo" + " para " + cov_value
                + "\u03C3" + " en Z"),
            color="orange")
        plt.plot(0, 0, "ko", label="Punto de partida")
        plt.plot(
            x_values_pose[frame_beg], z_values_pose[frame_beg],
            "co",
            label=("Posición estimada del vehículo al"
                   + " volver al inicio de trayectoria"
                   + " para \u03C3 original."))
        plt.plot(
            x_values_pose_1[frame_beg], z_values_pose_1[frame_beg],
            "co", color="purple",
            label=("Posición estimada del vehículo al"
                   + " volver al inicio de trayectoria"
                   + " para " + cov_value + "\u03C3."))
        plt.ylabel("Eje Z [m]")
        plt.xlabel("Eje X [m]")
        plt.xlim(-2.5, 0.5)
        plt.ylim(-3, 7)
        plt.legend(loc='upper left')
        plt.title("Trayectoria del vehículo")

    asdasd = (6, 5)
    if plot_dict["TrayectoriaX"]:
        plt.figure(figsize=asdasd)
        plt.plot(lista_idxs, x_values_GT, label="Ground truth en eje X")
        plt.plot(
            lista_idxs,
            x_values_odometry, "--",
            label=("Dead reckoning para " + "\u03C3" + " original en X"),
            color="green")
        plt.plot(
            lista_idxs,
            x_values_odometry_1,
            label=(
                "Dead reckoning" + " para " + cov_value + "\u03C3" + " en X"),
            color="darkgreen")
        plt.plot(
            lista_idxs,
            x_values_pose, "--",
            label=(
                "Pose estimada del vehículo para " + "\u03C3"
                + " original en X"),
            color="red")
        plt.plot(
            lista_idxs,
            x_values_pose_1,
            label=(
                "Pose estimada del vehículo" + " para " + cov_value
                + "\u03C3" + " en X"),
            color="orange")
        plt.plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
        plt.ylabel("Distancia del origen [m]")
        plt.xlabel("Tiempo [frames]")
        plt.legend(loc='upper left')
        plt.ylim(-2.5, 1.5)
        plt.title("Comparación de trayectoria en el eje X")

    if plot_dict["TrayectoriaZ"]:
        plt.figure(figsize=asdasd)
        plt.plot(lista_idxs, z_values_GT, label="Ground truth en eje Z")
        plt.plot(
            lista_idxs,
            z_values_odometry, "--",
            label=("Dead reckoning para " + "\u03C3" + " original en Z"),
            color="green")
        plt.plot(
            lista_idxs,
            z_values_odometry_1,
            label=(
                "Dead reckoning" + " para " + cov_value + "\u03C3" + " en Z"),
            color="darkgreen")
        plt.plot(
            lista_idxs,
            z_values_pose, "--",
            label=(
                "Pose estimada del vehículo para " + "\u03C3"
                + " original en Z"),
            color="red")
        plt.plot(
            lista_idxs,
            z_values_pose_1,
            label=(
                "Pose estimada del vehículo" + " para " + cov_value
                + "\u03C3" + " en Z"),
            color="orange")
        plt.plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
        plt.ylabel("Distancia del origen [m]")
        plt.xlabel("Tiempo [frames]")
        plt.legend(loc='upper left')
        plt.ylim(-2.5, 7.5)
        plt.title("Comparación de trayectoria en el eje Z")

    if plot_dict["TrayectoriaY"]:
        plt.figure(figsize=asdasd)
        plt.plot(
            lista_idxs, [0]*len(lista_idxs),
            label="Ground truth en eje Y")
        plt.plot(
            lista_idxs,
            df_odometry.y, "--",
            label=("Dead reckoning para " + "\u03C3" + " original en Y"),
            color="green")
        plt.plot(
            lista_idxs,
            df_odometry_1.y,
            label=(
                "Dead reckoning" + " para " + cov_value + "\u03C3" + " en Y"),
            color="darkgreen")
        plt.plot(
            lista_idxs,
            filtered_df_pose.y, "--",
            label=(
                "Pose estimada del vehículo para " + "\u03C3"
                + " original en Y"),
            color="red")
        plt.plot(
            lista_idxs,
            filtered_df_pose_1.y,
            label=(
                "Pose estimada del vehículo" + " para " + cov_value
                + "\u03C3" + " en Y"),
            color="orange")
        plt.plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
        plt.ylabel("Distancia del origen [m]")
        plt.xlabel("Tiempo [frames]")
        plt.legend(loc='best')
        plt.title("Comparación de trayectoria en el eje Y")
    plt.show()
