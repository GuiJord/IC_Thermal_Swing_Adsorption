import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator

Rg = 8.314e-3                               #kJ/mol.K
# T = 45 + 273.15                              #K
N = 100

partial_p_interval = np.linspace(0,1,N)     #bar

#Parameters
b0 = np.array([3.33e-6,3.21e-5,6.54e-13])   #1/bar^m
dH = -np.array([41.61,18.79,95.66])         #kJ/mol
dH_K = -np.array([21.25,13.81,46.53])       #kJ/mol
m = np.array([1.21,1.25,1.969])             #

#Pure Component Isotherms
def calc_loading_single(T,partial_p,calf_type):
    if calf_type == "crystal":
        #Crystals
        K0 = np.array([1.53e-4,6.96e-4,5.47e-7])    #mmol/g crystal
        q_s = np.array([2.959,2.959,7.603])         #mmol/g crystal
    elif calf_type == "particle":
        #Particles
        K0 = np.array([1.45e-4,6.58e-4,5.17e-7])    #mmol/g particle
        q_s = np.array([2.799,2.799,7.192])         #mmol/g particle

    b = b0*np.exp(-dH/(Rg*T))
    K = K0*np.exp(-dH_K/(Rg*T))
    q_star = q_s*b*partial_p**m/(1+b*partial_p**m) + K*partial_p
    return q_star


#Mixture Isotherms
# calf_type = "particle"
#Particles
# K0 = np.array([1.45e-4,6.58e-4,5.17e-7])
# q_s = np.array([2.799,2.799,7.192])

A = np.array([0.254,2.9e8])
B = np.array([3.49,0.044])
A_prime = 8.38e6
B_prime = 0.05
E = 97.7

def calc_loading_mixture(T,P,comp,RH,calf_type):    
    if calf_type == "crystal":
        #Crystals
        K0 = np.array([1.53e-4,6.96e-4,5.47e-7])    #mmol/g crystal
        q_s = np.array([2.959,2.959,7.603])         #mmol/g crystal
    elif calf_type == "particle":
        #Particles
        K0 = np.array([1.45e-4,6.58e-4,5.17e-7])    #mmol/g particle
        q_s = np.array([2.799,2.799,7.192])         #mmol/g particle


    partial_p = P*comp
    b = b0*np.exp(-dH/(Rg*T))
    K = K0*np.exp(-dH_K/(Rg*T))

    d_CO2 = 101.33
    d_N2 = np.exp(-E*partial_p[1]/T)
    d_H2O = A_prime*(partial_p[0])**2*np.exp(-B_prime*T)
    
    aux = (RH/A[0])**B[0]
    gamma_CO2 = 1 - (aux/(1+aux))**(1-RH**d_CO2)
    gamma_N2 = 1
    aux = (RH/A[1])**B[1]
    gamma_H2O = d_N2*(aux/(1+aux))**(1-RH**d_H2O)
    
    gamma = np.array([gamma_CO2,gamma_N2,gamma_H2O])
    
    q_star_mix = q_s*b*gamma*partial_p**m/(1+np.sum(b*gamma*partial_p**m)) + K*partial_p
    return q_star_mix

#Psat coeff 293 to 343 K
A_psat = 6.20963
B_psat = 2354.731
C_psat = 7.559
def calc_Psat(T):
    Psat = 10**(A_psat-B_psat/(T+C_psat))
    return Psat


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



# fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)
# P = 1 #bar
# comp_interval = np.linspace(0,1,N)
# RH_interval = np.linspace(0,1,N)
# T_interval = [40,60]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15
#     Psat = calc_Psat(T)
#     # y_H2O_inverval = RH_interval*Psat/P
#     q_star_array = np.zeros((3,N))
#     for i in range(0,N):
#         # y_CO2 = 1-y_H2O_inverval[i]
#         # y_N2 = 0
#         # y_H2O = y_H2O_inverval[i]
        
#         y_CO2 = comp_interval[i]
#         y_N2 = 1 - y_CO2
#         y_H2O = 0
        
#         comp = np.array([y_CO2,y_N2,y_H2O])
#         RH = 0
#         q_star_array[:,i] = calc_loading_mixture(T,P,comp,RH)
#         plt.plot(partial_p_interval,q_star_array[1],label=f'{round(T-273.15)}°C',color=tol[j])


#Fig 2a
# P = 1 #bar
# comp_interval = np.linspace(0,1,N)

# fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)
# ax2 = ax1.twinx()

# calf_type = 'particle'

# T_interval = [40]
# for j, T_C in enumerate(T_interval):
#     T = T_C + 273.15
#     Psat = calc_Psat(T)
#     q_star_array = np.zeros((3, N))

#     for i in range(N):
#         y_CO2 = comp_interval[i]
#         y_N2 = 1 - y_CO2
#         y_H2O = 0

#         comp = np.array([y_CO2, y_N2, y_H2O])
#         RH = 0
#         q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH, calf_type)

#     # CO2
#     ax1.plot(comp_interval,q_star_array[0],color=tol[j],linestyle="-",label=f"CO₂ (particle)")
#     # N2
#     ax2.plot(comp_interval,q_star_array[1],color=tol[j+1],linestyle="-",label=f"N₂ (particle)")

# calf_type = 'crystal'

# T_interval = [40]
# for j, T_C in enumerate(T_interval):
#     T = T_C + 273.15
#     Psat = calc_Psat(T)
#     q_star_array = np.zeros((3, N))

#     for i in range(N):
#         y_CO2 = comp_interval[i]
#         y_N2 = 1 - y_CO2
#         y_H2O = 0

#         comp = np.array([y_CO2, y_N2, y_H2O])
#         RH = 0
#         q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH, calf_type)

#     # CO2
#     ax1.plot(comp_interval,q_star_array[0],color=tol[j],linestyle="--",label=f"CO₂ (crystal)")
#     # N2
#     ax2.plot(comp_interval,q_star_array[1],color=tol[j+1],linestyle="--",label=f"N₂ (crystal)")

# article_CO2 = np.load('CO2_particle_40.npy')
# article_N2 = np.load('N2_particle_40.npy')

# ax1.plot(article_CO2[:,0],article_CO2[:,1],color='black',linestyle=':',ms=3,marker='s',label='artigo')
# ax2.plot(article_N2[:,0],article_N2[:,1],color='black',linestyle=':',ms=3,marker='s')

# ax1.set_box_aspect(1)

# ax1.set_xlabel("CO₂ mole fraction (-)")
# ax1.set_ylabel(f"CO₂ loading (mmol/g)")
# ax1.set_xlim(0, 1)
# ax1.set_ylim(0, 4)

# ax2.set_ylabel(f"N₂ loading (mmol/g)")
# ax2.set_ylim(0, 0.3)

# # Legenda conjunta
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.legend(lines1 + lines2, labels1 + labels2,loc='center right')

# plt.tight_layout()
# plt.savefig("isoterma_mix_Fig2a.png", dpi=200)



#Fig 2b
P = 1 #bar
comp_interval = np.linspace(0,1,N)

fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)
ax2 = ax1.twinx()

calf_type = 'particle'

T_interval = [60]
for j, T_C in enumerate(T_interval):
    T = T_C + 273.15
    Psat = calc_Psat(T)
    q_star_array = np.zeros((3, N))

    for i in range(N):
        y_CO2 = comp_interval[i]
        y_N2 = 1 - y_CO2
        y_H2O = 0

        comp = np.array([y_CO2, y_N2, y_H2O])
        RH = 0
        q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH, calf_type)

    # CO2
    ax1.plot(comp_interval,q_star_array[0],color=tol[j],linestyle="-",label=f"CO₂ (particle)")
    # N2
    ax2.plot(comp_interval,q_star_array[1],color=tol[j+1],linestyle="-",label=f"N₂ (particle)")

calf_type = 'crystal'

T_interval = [60]
for j, T_C in enumerate(T_interval):
    T = T_C + 273.15
    Psat = calc_Psat(T)
    q_star_array = np.zeros((3, N))

    for i in range(N):
        y_CO2 = comp_interval[i]
        y_N2 = 1 - y_CO2
        y_H2O = 0

        comp = np.array([y_CO2, y_N2, y_H2O])
        RH = 0
        q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH, calf_type)

    # CO2
    ax1.plot(comp_interval,q_star_array[0],color=tol[j],linestyle="--",label=f"CO₂ (crystal)")
    # N2
    ax2.plot(comp_interval,q_star_array[1],color=tol[j+1],linestyle="--",label=f"N₂ (crystal)")

article_CO2 = np.load('CO2_particle_60.npy')
article_N2 = np.load('N2_particle_60.npy')

ax1.plot(article_CO2[:,0],article_CO2[:,1],color='black',linestyle=':',ms=3,marker='s',label='artigo')
ax2.plot(article_N2[:,0],article_N2[:,1],color='black',linestyle=':',ms=3,marker='s')

ax1.set_box_aspect(1)

ax1.set_xlabel("CO₂ mole fraction (-)")
ax1.set_ylabel(f"CO₂ loading (mmol/g)")
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 3.5)

ax2.set_ylabel(f"N₂ loading (mmol/g)")
ax2.set_ylim(0, 0.2)

# Legenda conjunta
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax1.legend(lines1 + lines2, labels1 + labels2,loc='center right')

plt.tight_layout()
plt.savefig("isoterma_mix_Fig2b.png", dpi=200)



# #Fig 4a
# P = 0.97 #bar
# RH_interval = np.linspace(0,1,N)

# fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)

# ax2 = ax1.twinx()

# T_interval = [22]

# for j, T_C in enumerate(T_interval):

#     T = T_C + 273.15
#     Psat = calc_Psat(T)

#     q_star_array = np.zeros((3, N))

#     y_H2O_inverval = RH_interval*Psat/P
#     for i in range(N):
#         y_CO2 = 1-y_H2O_inverval[i]
#         y_N2 = 0
#         y_H2O = y_H2O_inverval[i]

#         comp = np.array([y_CO2, y_N2, y_H2O])

#         q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH_interval[i])

#     # CO2
#     ax1.plot(
#         RH_interval*100,
#         q_star_array[0],
#         color=tol[j],
#         linestyle="-",
#         label=f"CO₂ ({T_C}°C)"
#     )

#     # N2
#     ax2.plot(
#         RH_interval*100,
#         q_star_array[2],
#         color=tol[j+1],
#         linestyle="-",
#         label=f"H₂O ({T_C}°C)"
#     )

# ax1.set_box_aspect(1)

# ax1.set_xlabel("RH (%)")
# ax1.set_ylabel(f"CO₂ loading (mmol/g {calf_type})")
# ax1.set_xlim(0, 75)
# ax1.set_ylim(0, 8)

# ax2.set_ylabel(f"H₂O loading (mmol/g {calf_type})")
# ax2.set_ylim(0, 8)

# # Legenda conjunta
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.legend(lines1 + lines2, labels1 + labels2,loc="upper left")

# plt.tight_layout()
# plt.savefig(f"isoterma_mix_Fig4a_{calf_type}.png", dpi=200)



# #Fig 4b
# P = 0.15 #bar
# RH_interval = np.linspace(0,1,N)

# fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)

# ax2 = ax1.twinx()

# T_interval = [25,35,45]

# for j, T_C in enumerate(T_interval):

#     T = T_C + 273.15
#     Psat = calc_Psat(T)

#     q_star_array = np.zeros((3, N))

#     y_H2O_inverval = RH_interval*Psat/P
#     for i in range(N):
#         y_CO2 = 1-y_H2O_inverval[i]
#         y_N2 = 0
#         y_H2O = y_H2O_inverval[i]

#         comp = np.array([y_CO2, y_N2, y_H2O])

#         q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH_interval[i])

#     # CO2
#     ax1.plot(
#         RH_interval*100,
#         q_star_array[0],
#         color=tol[j],
#         linestyle="-",
#         label=f"CO₂ ({T_C}°C)"
#     )

#     # N2
#     ax2.plot(
#         RH_interval*100,
#         q_star_array[2],
#         color=tol[j],
#         linestyle="--",
#         label=f"H₂O ({T_C}°C)"
#     )

# ax1.set_box_aspect(1)

# ax1.set_xlabel("RH (%)")
# ax1.set_ylabel(f"CO₂ loading (mmol/g {calf_type})")
# ax1.set_xlim(0, 60)
# ax1.set_ylim(0, 8)

# ax2.set_ylabel(f"H₂O loading (mmol/g {calf_type})")
# ax2.set_ylim(0, 8)

# # Legenda conjunta
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.legend(lines1 + lines2, labels1 + labels2,loc="upper left")

# plt.tight_layout()
# plt.savefig(f"isoterma_mix_Fig4b_{calf_type}.png", dpi=200)



# #Fig 5b
# P = 1 #bar
# RH_interval = np.linspace(0,1,N)

# fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)

# ax2 = ax1.twinx()

# T_interval = [40]

# for j, T_C in enumerate(T_interval):

#     T = T_C + 273.15
#     Psat = calc_Psat(T)

#     q_star_array = np.zeros((3, N))

#     y_H2O_inverval = RH_interval*Psat/P
#     for i in range(N):
#         y_CO2 = 0
#         y_N2 = 1-y_H2O_inverval[i]
#         y_H2O = y_H2O_inverval[i]

#         comp = np.array([y_CO2, y_N2, y_H2O])

#         q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH_interval[i])

#     # H2O
#     ax1.plot(
#         RH_interval*100,
#         q_star_array[2],
#         color=tol[j+1],
#         linestyle="-",
#         label=f"H₂O ({T_C}°C)"
#     )

#     # N2
#     ax2.plot(
#         RH_interval*100,
#         q_star_array[1],
#         color=tol[j],
#         linestyle="-",
#         label=f"N₂ ({T_C}°C)"
#     )

# ax1.set_box_aspect(1)

# ax1.set_xlabel("RH (%)")
# ax1.set_ylabel(f"H₂O loading (mmol/g {calf_type})")
# ax1.set_xlim(0, 100)
# ax1.set_ylim(0, 10)

# ax2.set_ylabel(f"N₂ loading (mmol/g {calf_type})")
# ax2.set_ylim(0, 0.3)

# # Legenda conjunta
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.legend(lines1 + lines2, labels1 + labels2,loc="upper left")

# plt.tight_layout()
# plt.savefig(f"isoterma_mix_Fig5b_{calf_type}.png", dpi=200)




# #Fig 6b
# P = 1 #bar
# RH_interval = np.linspace(0,1,N)

# fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)

# ax2 = ax1.twinx()

# T_interval = [30]

# for j, T_C in enumerate(T_interval):

#     T = T_C + 273.15
#     Psat = calc_Psat(T)

#     q_star_array = np.zeros((3, N))

#     y_H2O_inverval = RH_interval*Psat/P
#     for i in range(N):
#         y_CO2 = 0.3
#         y_N2 = 1-y_CO2-y_H2O_inverval[i]
#         y_H2O = y_H2O_inverval[i]

#         comp = np.array([y_CO2, y_N2, y_H2O])

#         q_star_array[:, i] = calc_loading_mixture(T, P, comp, RH_interval[i])

#     # CO2
#     ax1.plot(
#         RH_interval*100,
#         q_star_array[0],
#         color=tol[j],
#         linestyle="-",
#         label=f"CO₂ ({T_C}°C)"
#     )

#     # H2O
#     ax1.plot(
#         RH_interval*100,
#         q_star_array[2],
#         color=tol[j+2],
#         linestyle="-",
#         label=f"H₂O ({T_C}°C)"
#     )

#     # N2
#     ax2.plot(
#         RH_interval*100,
#         q_star_array[1],
#         color=tol[j+1],
#         linestyle="-",
#         label=f"N₂ ({T_C}°C)"
#     )

# ax1.set_box_aspect(1)

# ax1.set_xlabel("RH (%)")
# ax1.set_ylabel(f"CO₂/H₂O loading (mmol/g {calf_type})")
# ax1.set_xlim(0, 100)
# ax1.set_ylim(0, 9)

# ax2.set_ylabel(f"N₂ loading (mmol/g {calf_type})")
# ax2.set_ylim(0, 0.2)
# # ax2.yaxis.set_major_locator(MultipleLocator(0.05))
# ax2.yaxis.set_major_locator(MaxNLocator(5))

# # Legenda conjunta
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# ax1.legend(lines1 + lines2, labels1 + labels2,loc="upper left")

# plt.tight_layout()
# plt.savefig(f"isoterma_mix_Fig6b_{calf_type}.png", dpi=200)