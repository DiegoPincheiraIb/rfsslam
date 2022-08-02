import pandas as pd
import numpy as np
import collections
import copy
import matplotlib.pyplot as plt


df_pose = pd.read_table(
    "/home/diego/RFS_SLAM/rfsslam/build/data/rbphdslam6d/particlePose.dat",
    sep="   ",
    names=[
        "timestamp", "particle",
        "x", "y", "z", "qx", "qy", "qz", "qw", "weight"])

df_odometry = pd.read_table(
    "/home/diego/RFS_SLAM/rfsslam/build/data/rbphdslam6d/deadReckoning.dat",
    names=["timestamp", "x", "y", "z", "qx", "qy", "qz", "qw"], sep="   ")

del df_pose["particle"]
del df_pose["qx"]
del df_pose["qz"]


del df_odometry["qx"]
del df_odometry["qz"]

df_pose["timestamp"] = df_pose["timestamp"].div(3)

# print("df pose: \n", df_pose.head(), "\n")
# print("df odometry: \n", df_odometry)
print_info = False
if print_info:
    print("lista de timestamps en df_pose",
          list(collections.Counter(df_pose["timestamp"])))

    print("lista de timestamps en df_odometry",
          list(collections.Counter(df_odometry["timestamp"])))

df_pose["weight"].fillna(1, inplace=True)

df_pose = df_pose.drop_duplicates()
lista_tstmp = list(collections.Counter(df_pose["timestamp"]))
asdasd = []
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
    asdasd.append(media_obj)
    # print(np.subtract(np.array(media_obj), np.array(df_odometry)))

del df_odometry["timestamp"]
# del df_odometry["timestamp"]

new_odometry = copy.deepcopy(df_odometry).cumsum()
new_odometry["qw"] = df_odometry["qw"]
# print(new_odometry.head())

filtered_df_pose = pd.DataFrame(
        asdasd, columns=["x", "y", "z", "qy", "qw"])
# print(filtered_df_pose)
# print(df_odometry)
error_df = df_odometry.subtract(filtered_df_pose)
# print(error_df["x"].isnull().sum())
# error_df.plot(y=["x", "z"], label=["Error en eje X", "Error en eje Z"])
lista_idxs = list(error_df.index)


x_values_odometry = list(df_odometry["x"])
y_values_odometry = list(df_odometry["y"])
z_values_odometry = list(df_odometry["z"])
x_values_pose = list(filtered_df_pose["x"])
y_values_pose = list(filtered_df_pose["y"])
z_values_pose = list(filtered_df_pose["z"])
x_values_odometry = [value * -1 for value in x_values_odometry]
x_values_pose = [value * -1 for value in x_values_pose]
# print(x_values_odometry)
# print(z_values_odometry)
# print(x_values_pose)
# print(z_values_pose)

frame_beg = 180

# plt.plot(lista_idxs, z_values_odometry)
# plt.plot(lista_idxs, z_values_pose)
# plt.plot(x_values_pose, z_values_pose)
# plt.plot(x_values_odometry, z_values_odometry)
hola = True
if hola:
    print(
        "error en ejes: ",
        error_df.loc[frame_beg, "x"],
        error_df.loc[frame_beg, "y"],
        error_df.loc[frame_beg, "z"],
        error_df.loc[frame_beg, "qy"],
        error_df.loc[frame_beg, "qw"])
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
asss = True
if asss:
    plt.plot(lista_idxs, error_df.x, label="Error en el eje X")
    plt.plot(lista_idxs, error_df.y, label="Error en el eje Y")
    plt.plot(lista_idxs, error_df.z, label="Error en el eje Z")
    plt.plot(frame_beg, 0, "go", label="Retorno a la posición inicial")
    plt.ylabel("Error [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.ylim(-0.7, 1.2)
    plt.title("Error estimado a lo largo de la trayectoria")

    plt.figure()
    plt.plot(
        x_values_odometry, z_values_odometry, label="Ground truth")
    plt.plot(x_values_pose, z_values_pose, label="Pose estimada del vehículo")
    # plt.plot(lista_idxs, z_values_odometry)
    # plt.plot(lista_idxs, z_values_pose)
    plt.plot(0, 0, "go", label="Punto de partida")
    plt.plot(
        x_values_pose[frame_beg], z_values_pose[frame_beg],
        "ro",
        label="Posición estimada del vehículo al volver al inicio de trayectoria")
    plt.ylabel("Eje Z [m]")
    plt.xlabel("Eje X [m]")
    plt.xlim(-3, 1)
    plt.ylim(-3, 7)
    plt.legend(loc='upper left')
    plt.title("Trayectoria del vehículo")

    plt.figure()
    plt.plot(lista_idxs, x_values_odometry, label="Ground truth en eje X")
    plt.plot(lista_idxs, x_values_pose, label="Pose estimada del vehículo en eje X")
    plt.plot(frame_beg, 0, "go", label="Retorno a la posición inicial")
    plt.ylabel("Distancia del origen [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.ylim(-2.5, 1.5)
    plt.title("Comparación de trayectoria en el eje X")

    plt.figure()
    plt.plot(lista_idxs, df_odometry.z, label="Ground truth en eje Z")
    plt.plot(
        lista_idxs, filtered_df_pose.z, label="Pose estimada del vehículo en eje Z")
    plt.plot(frame_beg, 0, "go", label="Retorno a la posición inicial")
    plt.ylabel("Distancia del origen [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.ylim(-2.5, 6.5)
    plt.title("Comparación de trayectoria en el eje Z")

    plt.figure()
    plt.plot(lista_idxs, [0]*len(lista_idxs), label="Ground truth en eje Y")
    plt.plot(
        lista_idxs, filtered_df_pose.y, label="Pose estimada del vehículo en eje Y")
    plt.plot(frame_beg, 0, "go", label="Retorno a la posición inicial")
    plt.ylabel("Distancia del origen [m]")
    plt.xlabel("Tiempo [frames]")
    plt.legend(loc='upper left')
    plt.title("Comparación de trayectoria en el eje Y")
    plt.show()
