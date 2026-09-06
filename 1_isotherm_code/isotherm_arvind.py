import numpy as np
import matplotlib.pyplot as plt

R = 8.314                                       #J/mol.K
# T = 45 + 273.15                               #K
N = 100

#Parameters
b0 = np.array([8.65*1e-7,2.5*1e-6])
d0 = np.array([2.63*1e-8,0])

dU_b = np.array([-36641.21,-1.58*1e4])
dU_d = np.array([-35690.66,0])

q_sb = np.array([3.09,5.84])
q_sd = np.array([2.54,0])

#Pure Component Isotherms
def q_star_pure_calc(concentration, temperature):
    b = b0*np.exp(-dU_b/(R*temperature))
    d = d0*np.exp(-dU_d/(R*temperature))

    A1 = q_sb*b/(1+(b*concentration))
    A2 = q_sd*d/(1+(d*concentration))

    q_star = concentration*(A1 + A2)
    return q_star

def q_star_mix_calc(concentration, temperature):
    b = b0*np.exp(-dU_b/(R*temperature))
    d = d0*np.exp(-dU_d/(R*temperature))

    A1 = q_sb*b/(1+np.sum(b*concentration))
    A2 = q_sd*d/(1+np.sum(d*concentration))

    q_star = concentration*(A1 + A2)
    return q_star


#Plotting
# cmap = plt.get_cmap("Paired")
tol = [
    "#4477AA",  # blue
    "#EE6677",  # red
    "#228833",  # green
    "#CCBB44",  # yellow
    "#66CCEE",  # cyan
    "#AA3377",  # purple
    "#BBBBBB",  # gray
]

#Fig 3
P_interval = np.linspace(0,0.15,N)
T = 25

plt.figure(figsize=(6, 4),dpi=200)
T_K = T + 273.15
q_star_array = np.zeros((2,N))
for i in range(0,N):
    concentration = P_interval[i]/(R*T_K)*1e5*np.array([1,0])
    q_star_array[:,i] = q_star_pure_calc(concentration, T_K)

plt.plot(P_interval,q_star_array[0],label=f'CO$_2$',color=tol[0])
for i in range(0,N):
    concentration = P_interval[i]/(R*T_K)*1e5*np.array([0,1])
    q_star_array[:,i] = q_star_pure_calc(concentration, T_K)

plt.plot(P_interval,q_star_array[1],label=f'N$_2$',color=tol[1])

q_CO2_arvind = np.load('CO2_arvind.npy')
q_N2_arvind = np.load('N2_arvind.npy')
plt.plot(q_CO2_arvind[:,0],q_CO2_arvind[:,1],marker='s',ms=3,ls='none',color='black', label='artigo')
plt.plot(q_N2_arvind[:,0],q_N2_arvind[:,1],marker='s',ms=3,ls='none',color='black')


plt.xlabel("Pressure (bar)")
plt.ylabel(f"Loading (mmol/g)")
plt.xlim(0,0.15)
plt.ylim(0,4)
plt.legend()
plt.tight_layout()
plt.savefig('isoterma_Fig3_Arvind_pure.png',dpi=200)

#Fig 
P_interval = np.linspace(0,0.15,N)
T = 25

plt.figure(figsize=(6, 4),dpi=200)
T_K = T + 273.15
q_star_array = np.zeros((2,N))
for i in range(0,N):
    concentration = P_interval[i]/(R*T_K)*1e5*np.array([0.15,0.85])
    q_star_array[:,i] = q_star_mix_calc(concentration, T_K)

plt.plot(P_interval,q_star_array[0],label=f'CO$_2$ mix',color=tol[0])
plt.plot(P_interval,q_star_array[1],label=f'N$_2$ mix',color=tol[1])

q_star_array = np.zeros((2,N))
for i in range(0,N):
    concentration = P_interval[i]/(R*T_K)*1e5*np.array([1,0])
    q_star_array[:,i] = q_star_pure_calc(concentration, T_K)

plt.plot(P_interval,q_star_array[0],ls='--',label=f'CO$_2$ pure',color=tol[0])
for i in range(0,N):
    concentration = P_interval[i]/(R*T_K)*1e5*np.array([0,1])
    q_star_array[:,i] = q_star_pure_calc(concentration, T_K)

plt.plot(P_interval,q_star_array[1],ls='--',label=f'N$_2$ pure',color=tol[1])

plt.xlabel("Pressure (bar)")
plt.ylabel(f"Loading (mmol/g)")
plt.xlim(0,0.15)
plt.ylim(0,4)
plt.legend()
plt.tight_layout()
plt.savefig('isoterma_Fig_Arvind_mix.png',dpi=200)
