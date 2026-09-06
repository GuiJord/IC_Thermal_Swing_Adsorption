import numpy as np
import pandas as pd

nome = str(input("Digite o nome do .CSV: "))

# Lê o CSV do WebPlotDigitizer
df = pd.read_csv(
    f"{nome}.csv",
    sep=";",
    decimal=","
)

# Converte para um numpy.array
data = df.to_numpy()

# Salva em formato .npy
np.save(f"{nome}.npy", data)