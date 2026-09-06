import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter, NullFormatter, FixedLocator


N_list = np.array([10,20,30,40,50,60,70,80,90,100,200,300,400,500])

time_list = np.array([])

for N in N_list:
    with open(f"./run_{N}/results.txt") as f:
        numbers = f.read().split()
        time = float(numbers[3])
    time_list = np.append(time_list,time)

time_list = time_list/min(time_list)

fig, ax1 = plt.subplots(figsize=(9,6), dpi=200)

# Error MB
ax1.plot(
    N_list,
    time_list,
    color='green',
    linestyle='-',
    ms=6,marker='o',
    label='t$_i$/t$_{10}$'
)

ax1.set_box_aspect(6/9)

# Labels
ax1.set_xlabel("Número de Volumes Finitos [-]", fontsize=13)
ax1.set_ylabel("Tempo Computacional [-]", fontsize=13)

# ax1.set_xscale("log")
ax1.set_yscale("log")


    # Tick appearance
ax1.tick_params(axis='both', which='major',
                direction='in', length=7, width=1.2,
                labelsize=12, top=True)

ax1.tick_params(axis='both', which='minor',
                direction='in', length=4, width=0.8,
                top=True)

# Grid
ax1.grid(which='major', linestyle='--', alpha=0.3)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
ax1.legend(lines1, labels1, fontsize=12,loc='lower right')

plt.tight_layout()
plt.savefig("time_graph.png", dpi=200)