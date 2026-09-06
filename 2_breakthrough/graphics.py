import numpy as np
import matplotlib.pyplot as plt
from parameters import *

def plot_P_time():
    plt.figure(figsize=(9, 4))
    plt.plot(t, P_ad[-1,:],color='black')
    plt.ylabel('Pressure [-]')
    plt.xlabel('Time [-]')
    plt.xlim(left=0)
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('pressure_time.png',dpi=200)
    plt.close()

def plot_P_column():
    plt.figure(figsize=(9, 4))
    plt.plot(z, P_ad[:,-1],color='black')
    plt.ylabel('Pressure [-]')
    plt.xlabel('Z [-]')
    plt.xlim(left=0)
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('pressure_column.png',dpi=200)
    plt.close()

def plot_Tw_time():
    plt.figure(figsize=(9, 4))
    plt.plot(t, Tw_ad[-1,:],color='lime')
    plt.ylabel('Tw [-]')
    plt.xlabel('Time [-]')
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig('temperature_wall_time.png',dpi=200)
    plt.close()

def plot_Tw_column():
    plt.figure(figsize=(9, 4))
    plt.plot(z, Tw_ad[:,-1],color='lime')
    plt.ylabel('Tw [-]')
    plt.xlabel('Z [-]')
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig('temperature_wall_column.png',dpi=200)
    plt.close()

def plot_T_time():
    plt.figure(figsize=(8, 4))
    plt.plot(t, T_ad[-1,:],color='lime')
    plt.plot(T_arvind[:,0],T_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
    plt.ylabel('T [-]')
    plt.xlabel('Time [-]')
    plt.ylim(0.8, 1.4)
    plt.xlim(0, 8000)
    plt.xlim(left=0)
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('temperature_time.png',dpi=200)
    plt.close()

def plot_T_column():
    plt.figure(figsize=(9, 4))
    plt.plot(z, T_ad[:,-1],color='lime')
    plt.ylabel('T [-]')
    plt.xlabel('Z [-]')
    # plt.ylim(0.96, 1.04)
    # plt.xlim(0,1)
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('temperature_column.png',dpi=200)
    plt.close()

def plot_v_time():
    plt.figure(figsize=(9, 4))
    plt.plot(t, v_ad[:],color='dodgerblue')
    plt.plot(v_arvind[:,0],v_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
    plt.ylabel('v [-]')
    plt.xlabel('Time [-]')
    plt.ylim(0.8, 1.4)
    plt.xlim(0, 8000)
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig('velocity_time.png',dpi=200)
    plt.close()

def plot_v_column():
    plt.figure(figsize=(9, 4))
    # plt.plot(z, v_ad[:,-1],color='dodgerblue')
    plt.ylabel('v [-]')
    plt.xlabel('Z [-]')
    # plt.ylim(0.8, 1.4)
    # #plt.xlim(0.1, 15)
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig('velocity_column.png',dpi=200)
    plt.close()

def plot_x_time():
    plt.figure(figsize=(9, 4))
    plt.plot(t, x_CO2[-1,:],color='r', label='CO2')
    plt.plot(t, x_N2[-1,:],color='turquoise', label='N2')
    plt.ylabel('x [-]')
    plt.xlabel('Time [-]')
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('x_time.png',dpi=200)
    plt.close()

def plot_x_column():
    plt.figure(figsize=(9, 4))
    plt.plot(z, x_CO2[:,-1],color='r', label='CO2')
    plt.plot(z, x_N2[:,-1],color='turquoise', label='N2')
    plt.ylabel('x [-]')
    plt.xlabel('Z [-]')
    plt.xlim(left=0)
    plt.legend()
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('x_column.png',dpi=200)
    plt.close()

def plot_y_time():
    y_N2 = 1 - y_CO2

    plt.figure(figsize=(9, 4))
    plt.plot(t, y_CO2[-1,:],color='r',lw=2,label='código próprio')
    plt.plot(y_CO2_arvind[:,0],y_CO2_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
    # plt.plot(t, y_N2[-1,:],color='turquoise', label='N2')
    plt.ylabel('y$_{CO2}$')
    plt.xlabel('Time [-]')
    plt.xlim(0,8000)
    plt.ylim(0,0.16)
    # plt.locator_params(axis='x', nbins=)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('y_time.png',dpi=200)
    plt.close()

    
    plt.figure(figsize=(5, 4))
    plt.plot(t, y_CO2[-1,:],color='r',lw=2)#,label='CO2')
    plt.plot(y_peak1_arvind[:,0],y_peak1_arvind[:,1],marker='s',ms=5,ls='none',color='black', label='artigo')
    # plt.plot(t, y_N2[-1,:],color='turquoise', label='N2')
    plt.ylabel('y$_{CO2}$')
    plt.xlabel('Time [-]')
    plt.xlim(490,550)
    plt.ylim(0,0.14)
    # plt.locator_params(axis='x', nbins=)
    plt.legend()
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('y_time_first_peak.png',dpi=200)
    plt.close()


    plt.figure(figsize=(5, 4))
    plt.plot(t, y_CO2[-1,:],color='r',lw=2)#,label='CO2')
    plt.plot(y_peak2_arvind[:,0],y_peak2_arvind[:,1],marker='s',ms=5,ls='none',color='black', label='artigo')
    # plt.plot(t, y_N2[-1,:],color='turquoise', label='N2')
    plt.ylabel('y$_{CO2}$')
    plt.xlabel('Time [-]')
    plt.xlim(3000,4400)
    plt.ylim(0.12,0.16)
    plt.locator_params(axis='y', nbins=5)
    # plt.locator_params(axis='x', nbins=)
    plt.legend()
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('y_time_second_peak.png',dpi=200)
    plt.close()


def plot_y_column():
    y_N2 = 1 - y_CO2

    plt.figure(figsize=(9, 4))
    plt.plot(z, y_CO2[:,-1],color='r', label='CO2')
    # plt.plot(z, y_N2[:,-1],color='turquoise', label='N2')
    plt.ylabel('y$_{CO2}$')
    plt.xlabel('Z [-]')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.legend()
    plt.tight_layout()
    plt.ticklabel_format(style='plain', useOffset=False)
    plt.savefig('y_column.png',dpi=200)
    plt.close()


# test_number = int(input('Test number: '))
# nth_cycle = int(input('Cycle: '))

# path = f'./testset/test_{test_number}/cycle_{nth_cycle}/'

# sufixes = ['press','ads','blow','evac']



def graph():
    # plot_P_time()
    # plot_P_column()
    # plot_Tw_time()
    # plot_Tw_column()
    plot_T_time()
    # plot_T_column()
    plot_v_time()
    # plot_v_column()
    plot_x_time()
    # plot_x_column()
    plot_y_time()
    # plot_y_column()


t = np.load('t.npy')
P_ad = np.load('P_ad.npy')
T_ad = np.load('T_ad.npy')
v_ad = np.load('v_ad.npy')
y_CO2 = np.load('y_CO2.npy')
x_CO2 = np.load('x_CO2.npy')
x_N2 = np.load('x_N2.npy')
Tw_ad = np.load('Tw_ad.npy')

N_time_steps = len(t)
N = len(P_ad[:,0])
z = np.linspace(0,1,N)

# y_CO2_arvind = np.load('y_arvind.npy')
# T_arvind = np.load('T_arvind.npy')
# v_arvind = np.load('v_arvind.npy')

# y_peak1_arvind = np.load('y_peak1_arvind.npy')
# y_peak2_arvind = np.load('y_peak2_arvind.npy')

# graph()



# fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharex=True)

# axs[0].plot(t, v_ad[:],color='dodgerblue')
# axs[0].plot(v_arvind[:,0],v_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
# axs[0].set_ylabel('v/v$_0$ [-]')
# axs[0].set_ylim(0.8, 1.4)
# axs[0].set_xlim(0, 8000)

# axs[1].plot(t, T_ad[-1,:],color='lime')
# axs[1].plot(T_arvind[:,0],T_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
# axs[1].set_ylabel('T/T$_0$ [-]')
# axs[1].set_ylim(0.8, 1.4)


# axs[2].plot(t, y_CO2[-1,:],color='r',lw=2)
# axs[2].plot(y_CO2_arvind[:,0],y_CO2_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
# axs[2].set_ylabel(r'$y_{CO_2}$ [-]')
# axs[2].set_ylim(0, 0.16)

# axs[3].plot(t, x_CO2[-1,:],color='turquoise', label='CO2')
# axs[3].plot(t, x_N2[-1,:],color='limegreen', label='N2')
# axs[3].plot([],[],color='black', label='artigo')
# axs[3].set_ylabel('q/q$_0$ [-]')
# axs[3].set_ylim(0, 0.65)
# axs[3].legend(loc='lower right')
# axs[3].set_xlabel('Time [-]')

# # ==========================================================
# # Formatting
# # ==========================================================
# for ax in axs:
#     ax.ticklabel_format(style='plain', useOffset=False)

# plt.tight_layout()
# plt.savefig('time_breakthrough.png', dpi=300)



start = 10
end = 400
every = 50

start -= 1
end += 1

time_list = [i for i in range(start,end,every)]
t_index_list = []




time_to_save = start
for i in range(start,N_time_steps):
    if round(t[i]) > time_to_save:
        t_index_list.append(i)
        time_to_save += every
    if round(t[i]) > end:
        break


fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# ==========================================================
# Temperature
# ==========================================================
colors = plt.cm.plasma(np.linspace(0.1, 0.5, len(t_index_list)))

for index, color in zip(t_index_list, colors):
    axs[0].plot(
        z, T_ad[:, index],
        color=color,
        lw=3,
        label=f't = {int(round(t[index],0))}'
    )

axs[0].set_ylabel('T/T$_0$ [-]')
axs[0].set_ylim(0.98, 1.22)
axs[0].set_xlim(0, 1)
axs[0].legend(loc='center left', bbox_to_anchor=(1.02, 0.5))

# ==========================================================
# Solid loading
# ==========================================================
colors1 = plt.cm.Blues(np.linspace(0.1, 0.5, len(t_index_list)))
colors2 = plt.cm.Greens(np.linspace(0.1, 0.5, len(t_index_list)))

for index, c1, c2 in zip(t_index_list, colors1, colors2):
    axs[1].plot(z, q_s0*x_CO2[:, index], color=c1, lw=3,label=f't = {int(round(t[index],0))}')
    axs[1].plot(z, q_s0*x_N2[:, index], color=c2, lw=3)

# Species legend
axs[1].plot([], [], color='navy', label='CO$_2$')
axs[1].plot([], [], color='limegreen', label='N$_2$')

axs[1].set_ylabel('q/q$_0$ [-]')
axs[1].set_ylim(0, 3.7)
axs[1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5))


# ==========================================================
# Gas composition
# ==========================================================
colors = plt.cm.Reds(np.linspace(0.1, 0.5, len(t_index_list)))

for index, color in zip(t_index_list, colors):
    axs[2].plot(
        z, y_CO2[:, index],
        color=color,
        lw=3,
        label=f't = {int(round(t[index],0))}'
    )

axs[2].set_ylabel(r'$y_{CO_2}$ [-]')
axs[2].set_ylim(0, 0.16)
axs[2].legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
axs[2].set_xlabel('Z [-]')


# ==========================================================
# Formatting
# ==========================================================
for ax in axs:
    ax.ticklabel_format(style='plain', useOffset=False)

plt.tight_layout()
plt.savefig('column_profiles1.png', dpi=200)



start = 550
end = 2000
every = 150

start -= 1
end += 1

time_list = [i for i in range(start,end,every)]
t_index_list = []


time_to_save = start
for i in range(start,N_time_steps):
    if round(t[i]) > time_to_save:
        t_index_list.append(i)
        time_to_save += every
    if round(t[i]) > end:
        break



fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# ==========================================================
# Temperature
# ==========================================================
colors = plt.cm.plasma(np.linspace(0.5, 1, len(t_index_list)))

for index, color in zip(t_index_list, colors):
    axs[0].plot(
        z, T_ad[:, index],
        color=color,
        lw=3,
        label=f't = {int(round(t[index],0))}'
    )

axs[0].set_ylabel('T/T$_0$ [-]')
axs[0].set_ylim(0.98, 1.22)
axs[0].set_xlim(0, 1)
axs[0].legend(loc='center left', bbox_to_anchor=(1.02, 0.5))

# ==========================================================
# Solid loading
# ==========================================================
colors1 = plt.cm.Blues(np.linspace(0.5, 1, len(t_index_list)))
colors2 = plt.cm.Greens(np.linspace(0.5, 1, len(t_index_list)))

for index, c1, c2 in zip(t_index_list, colors1, colors2):
    axs[1].plot(z, q_s0*x_CO2[:, index], color=c1, lw=3,label=f't = {int(round(t[index],0))}')
    axs[1].plot(z, q_s0*x_N2[:, index], color=c2, lw=3)

# Species legend
axs[1].plot([], [], color='navy', label='CO$_2$')
axs[1].plot([], [], color='limegreen', label='N$_2$')

axs[1].set_ylabel('q/q$_0$ [-]')
axs[1].set_ylim(0, 3.7)
axs[1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5))


# ==========================================================
# Gas composition
# ==========================================================
colors = plt.cm.Reds(np.linspace(0.5, 1, len(t_index_list)))

for index, color in zip(t_index_list, colors):
    axs[2].plot(
        z, y_CO2[:, index],
        color=color,
        lw=3,
        label=f't = {int(round(t[index],0))}'
    )

axs[2].set_ylabel(r'$y_{CO_2}$ [-]')
axs[2].set_ylim(0, 0.16)
axs[2].legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
axs[2].set_xlabel('Z [-]')


# ==========================================================
# Formatting
# ==========================================================
for ax in axs:
    ax.ticklabel_format(style='plain', useOffset=False)

plt.tight_layout()
plt.savefig('column_profiles2.png', dpi=200)


# plt.figure(figsize=(10, 4))
# colors = plt.cm.plasma(np.linspace(0.1, 1, len(t_index_list)))   # or plasma, turbo, inferno, etc.

# for index, color in zip(t_index_list, colors):
#     plt.plot(z, T_ad[:,index], color=color,label=f't = {int(round(t[index],0))}',marker='none',ls='-',lw=5)

# plt.ylabel('T [-]')
# plt.xlabel('Z [-]')
# # plt.ylim(bottom=0.98)
# plt.ylim(0.98, 1.25)
# plt.xlim(0,1)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.tight_layout()
# plt.ticklabel_format(style='plain', useOffset=False)
# plt.savefig('temperature_column_nova.png',dpi=200)
# plt.close()



# plt.figure(figsize=(10, 4))

# colors = plt.cm.Reds(np.linspace(0.1, 1, len(t_index_list)))   # or plasma, turbo, inferno, etc.

# for index, color in zip(t_index_list, colors):
#     plt.plot(z, y_CO2[:,index], color=color,label=f't = {int(round(t[index],0))}',marker='none',ls='-',lw=5)

# plt.ylabel('y CO2 [-]')
# plt.xlabel('Z [-]')
# # plt.ylim(bottom=0.98)
# plt.ylim(0, 0.15)
# plt.xlim(0,1)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.tight_layout()
# plt.ticklabel_format(style='plain', useOffset=False)
# plt.savefig('y_CO2_column_nova.png',dpi=200)
# plt.close()

# plt.figure(figsize=(10, 4))

# colors1 = plt.cm.Blues(np.linspace(0.1, 1, len(t_index_list)))   # or plasma, turbo, inferno, etc.
# colors2 = plt.cm.OrRd(np.linspace(0.1, 1, len(t_index_list)))   # or plasma, turbo, inferno, etc.

# for index, color1, color2 in zip(t_index_list, colors1, colors2):
#     plt.plot(z, q_s0*x_CO2[:,index], color=color1,label=f't = {int(round(t[index],0))}',marker='none',ls='-',lw=5)
#     plt.plot(z, q_s0*x_N2[:,index], color=color2,marker='none',ls='-',lw=5)

# plt.plot(100,100,color='blue',label='CO2')
# plt.plot(100,100,color='red',label='N2')

# plt.ylabel('q [mmol/g]')
# plt.xlabel('Z [-]')
# # plt.ylim(bottom=0.98)
# plt.ylim(0, 3.7)
# plt.xlim(0,1)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.tight_layout()
# plt.ticklabel_format(style='plain', useOffset=False)
# plt.savefig('q_CO2_column_nova.png',dpi=200)
# plt.close()