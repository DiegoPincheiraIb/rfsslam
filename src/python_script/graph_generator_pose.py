import pandas as pd
import numpy as np
import collections
import copy
import matplotlib.pyplot as plt

# Read files
df_pose = pd.read_table(
    "/home/diego/RFS_SLAM/rfsslam/build/data/rbphdslam6d/particlePose.dat",
    sep="   ",
    names=[
        "timestamp", "particle",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])

df_GT = pd.read_table(
    "/home/diego/RFS_SLAM/rfsslam/build/data/rbphdslam6d/gtPose.dat",
    names=["timestamp", "x", "y", "z", "qx", "qy", "qz", "qw"], sep="   ")

df_odometry = pd.read_table(
    "/home/diego/RFS_SLAM/rfsslam/build/data/rbphdslam6d/deadReckoning.dat",
    names=["timestamp", "x", "y", "z", "qx", "qy", "qz", "qw"], sep="   ")

# Delete unnecesary columns
del df_pose["particle"]
del df_pose["qx"]
del df_pose["qz"]

del df_GT["qx"]
del df_GT["qz"]

del df_odometry["qx"]
del df_odometry["qz"]

# Adjust timestamps
df_pose["timestamp"] = df_pose["timestamp"].div(3)


# Se imprime lista de timestamps
print_info = False
if print_info:
    print("lista de timestamps en df_pose",
          list(collections.Counter(df_pose["timestamp"])))

    print("lista de timestamps en df_odometry",
          list(collections.Counter(df_odometry["timestamp"])))

# Se llenan los pesos que tienen N/A
df_pose["weight"].fillna(1, inplace=True)

# Se limpia el dataframe para capturar la pose de la partícula con la mayor
# probabilidad.
df_pose = df_pose.drop_duplicates()
lista_tstmp = list(collections.Counter(df_pose["timestamp"]))
pose_filtered_data = []
for timestamp in lista_tstmp:
    timestamp = int(timestamp)
    tstmp3 = df_pose.loc[df_pose.loc[:, "timestamp"] == lista_tstmp[timestamp]]
    lolazo = tstmp3[tstmp3.weight == tstmp3.weight.max()]
    # print(tstmp3)
    media_obj = list(lolazo.mean())[1:-1]
    if np.isnan(np.array(media_obj)).any():
        print("timestamp: ", timestamp)
        print("aqui: ", list(df_pose.iloc[timestamp, :]))
    # print(media_obj)
    pose_filtered_data.append(media_obj)
    # print(np.subtract(np.array(media_obj), np.array(df_odometry)))

# Se descartan los timestamps
del df_odometry["timestamp"]
del df_GT["timestamp"]

new_GT = copy.deepcopy(df_GT).cumsum()
new_GT["qw"] = df_GT["qw"]

new_odometry = copy.deepcopy(df_odometry).cumsum()
new_odometry["qw"] = df_odometry["qw"]

filtered_df_pose = pd.DataFrame(
        pose_filtered_data, columns=["x", "y", "z", "qy", "qw"])

error_df_odometry = df_odometry.subtract(filtered_df_pose)
error_df_GT = df_GT.subtract(filtered_df_pose)
lista_idxs = list(error_df_odometry.index)

# Print errores
print_info_filtered_df_pose = False
if print_info_filtered_df_pose:
    print(filtered_df_pose)
    print(df_odometry)
    print(error_df_odometry["x"].isnull().sum())

# Convertir valores de odometría a listas para poder graficarlos
x_values_odometry = list(df_odometry["x"])
y_values_odometry = list(df_odometry["y"])
z_values_odometry = list(df_odometry["z"])

x_values_GT = list(df_GT["x"])
y_values_GT = list(df_GT["y"])
z_values_GT = list(df_GT["z"])

x_values_pose = list(filtered_df_pose["x"])
y_values_pose = list(filtered_df_pose["y"])
z_values_pose = list(filtered_df_pose["z"])
x_values_odometry = [value * -1 for value in x_values_odometry]
x_values_GT = [value * -1 for value in x_values_GT]
x_values_pose = [value * -1 for value in x_values_pose]

# Frame donde se retorna a la posición inicial
frame_beg = 180

print_console_error = True
if print_console_error:
    print(
        "error en ejes: ",
        error_df_odometry.loc[frame_beg, "x"],
        error_df_odometry.loc[frame_beg, "y"],
        error_df_odometry.loc[frame_beg, "z"],
        error_df_odometry.loc[frame_beg, "qy"],
        error_df_odometry.loc[frame_beg, "qw"])
    print(
        "Coord en ejes: ",
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

# Plot graphics
plot_graphs = True
if plot_graphs:
    fig_obj, ax_obj = plt.subplots(3, 1, sharex='col', sharey='row')
    ax_obj[0].plot(
        lista_idxs, [0]*len(lista_idxs), label="Ground truth")
    ax_obj[0].plot(lista_idxs, error_df_GT.x, label="Error de dead reckoning", color="green")
    ax_obj[0].plot(lista_idxs, error_df_odometry.x, label="Error de trayectoria estimada", color="red")
    ax_obj[0].plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
    ax_obj[0].set_ylabel("Error en eje X[m]")
    ax_obj[0].legend(loc='upper left')
    ax_obj[0].set_title("Error estimado a lo largo de la trayectoria")
    # plt.subplot(312)
    ax_obj[1].plot(
        lista_idxs, [0]*len(lista_idxs), label="Ground truth")
    ax_obj[1].plot(lista_idxs, error_df_GT.y, label="Error de dead reckoning", color="green")
    ax_obj[1].plot(lista_idxs, error_df_odometry.y, label="Error de trayectoria estimada", color="red")
    ax_obj[1].set_ylabel("Error en eje Y[m]")
    ax_obj[1].plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
    ax_obj[1].legend(loc='upper left')
    # plt.subplot(313)
    ax_obj[2].plot(
        lista_idxs, [0]*len(lista_idxs), label="Ground truth")
    ax_obj[2].plot(lista_idxs, error_df_GT.z, label="Error de dead reckoning", color="green")
    ax_obj[2].plot(lista_idxs, error_df_odometry.z, label="Error de trayectoria estimada", color="red")
    plt.subplots_adjust(hspace=.0)
    ax_obj[2].plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
    ax_obj[2].legend(loc='upper left')
    plt.ylabel("Error en eje Z[m]")
    plt.xlabel("Tiempo [frames]")

    plt.figure()
    plt.plot(
        x_values_GT, z_values_GT, label="Ground truth")
    plt.plot(
        x_values_odometry, z_values_odometry,
        label="Dead reckoning", color="green")
    plt.plot(
        x_values_pose, z_values_pose,
        label="Pose estimada del vehículo", color="red")
    plt.plot(0, 0, "ko", label="Punto de partida")
    plt.plot(
        x_values_pose[frame_beg], z_values_pose[frame_beg],
        "co",
        label=("Posición estimada del vehículo al"
               + "volver al inicio de trayectoria"))
    plt.ylabel("Eje Z [m]")
    plt.xlabel("Eje X [m]")
    plt.xlim(-3, 1)
    plt.ylim(-3, 7)
    plt.legend(loc='upper left')
    plt.title("Trayectoria del vehículo")

    plt.figure()
    plt.plot(lista_idxs, x_values_GT, label="Ground truth en eje X")
    plt.plot(
        lista_idxs,
        x_values_odometry,
        label="Dead reckoning en eje X", color="green")
    plt.plot(
        lista_idxs,
        x_values_pose,
        label=("Pose estimada del vehículo en eje X"), color="red")
    plt.plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
    plt.ylabel("Distancia del origen [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.ylim(-2.5, 1.5)
    plt.title("Comparación de trayectoria en el eje X")

    plt.figure()
    plt.plot(
        lista_idxs, df_GT.z,
        label="Ground truth en eje Z")
    plt.plot(
        lista_idxs, df_odometry.z,
        label="Dead reckoning en eje Z", color="green")
    plt.plot(
        lista_idxs,
        filtered_df_pose.z,
        label="Pose estimada del vehículo en eje Z", color="red")
    plt.plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
    plt.ylabel("Distancia del origen [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.ylim(-2.5, 6.5)
    plt.title("Comparación de trayectoria en el eje Z")

    plt.figure()
    plt.plot(lista_idxs, [0]*len(lista_idxs), label="Ground truth en eje Y")
    plt.plot(
        lista_idxs,
        df_odometry.y,
        label="Dead reckoning del vehículo en eje Y")
    plt.plot(
        lista_idxs,
        filtered_df_pose.y,
        label="Pose estimada del vehículo en eje Y")
    plt.plot(frame_beg, 0, "ko", label="Retorno a la posición inicial")
    plt.ylabel("Distancia del origen [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.title("Comparación de trayectoria en el eje Y")
    plt.show()
