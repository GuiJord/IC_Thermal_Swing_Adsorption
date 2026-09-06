import numpy as np
import matplotlib.pyplot as plt

CSS_indicator_list = np.load('CSS_indicator_list.npy')
N = len(CSS_indicator_list)
cycles = np.linspace(1,N,N)

plt.plot(cycles,CSS_indicator_list,color='r')
plt.xlabel('Cycle')
plt.ylabel('CSS Indicator [-]')
plt.tight_layout()
plt.savefig('CSS_indicator.png',dpi=200)