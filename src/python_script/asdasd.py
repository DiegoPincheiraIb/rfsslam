import pandas as pd


df_obj = pd.DataFrame(columns=["dx", "dy", "dz", "dqx", "dqy", "dqz", "dqw"])
df_obj.loc[0] = [float(0) for _ in range(7)]
n_frames = 3
for frame in range(n_frames):
    valor_obj = []
    while len(valor_obj) != 7:
        valor_obj = input("Ingrese valores para x, y, z, i, j, k, w")
        valor_obj = valor_obj.split()
        valor_obj = [float(item_obj) for item_obj in valor_obj]
        if len(valor_obj) != 7:
            print("No se tienen 7 valores. Intente nuevamente")
    df_obj.loc[len(df_obj)] = valor_obj
print(df_obj)
