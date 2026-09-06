import numpy as np
import matplotlib.pyplot as plt

CSS_indicator_list = np.load('CSS_indicator_list.npy')
purity_list = np.load('purity_list.npy')
recovery_list = np.load('recovery_list.npy')

N = len(CSS_indicator_list)
cycles = np.linspace(1,N,N)

plt.figure(figsize=(6,5),dpi=200)
plt.plot(cycles,CSS_indicator_list,color='r')
plt.xlabel('Cycles [-]')
plt.ylabel('CSS Indicator [-]')
plt.tight_layout()
plt.savefig('CSS_indicator.png',dpi=200)
plt.close()

fig, ax1 = plt.subplots(figsize=(6,5),dpi=200)
ax2 = ax1.twinx()

# ax1.set_box_aspect(1)
ax1.plot(cycles,purity_list*100,color='royalblue',linestyle="--",label=f"Purity")
ax2.plot(cycles,recovery_list*100,color='limegreen',linestyle="-",label=f"Recovery")


ax1.set_xlabel("Cycles [-]")
ax1.set_ylabel(f"Purity [%]")
ax1.set_ylim(bottom=0,top=10)

ax2.set_ylabel(f"Recovery [%]")
ax2.set_ylim(bottom=0,top=40)

# Legenda conjunta
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2)

plt.tight_layout()
plt.savefig("purity_recovery_evolution.png", dpi=200)