import numpy as np
import matplotlib.pyplot as plt

Rg = 8.314e-3                               #kJ/mol.K
# T = 45 + 273.15                              #K
N = 100

partial_p_interval = np.linspace(0,1,N)     #bar

calf_type = "particle"# "crystal" #

#Parameters
b0 = np.array([3.33e-6,3.21e-5,6.54e-13])   #1/bar^m
dH = -np.array([41.61,18.79,95.66])         #kJ/mol
dH_K = -np.array([21.25,13.81,46.53])       #kJ/mol
m = np.array([1.21,1.25,1.969])             #

if calf_type == "crystal":
    #Crystals
    K0 = np.array([1.53e-4,6.96e-4,5.47e-7])    #mmol/g crystal
    q_s = np.array([2.959,2.959,7.603])         #mmol/g crystal
elif calf_type == "particle":
    #Particles
    K0 = np.array([1.45e-4,6.58e-4,5.17e-7])    #mmol/g particle
    q_s = np.array([2.799,2.799,7.192])         #mmol/g particle

#Pure Component Isotherms
def calc_loading_single(T,partial_p):
    b = b0*np.exp(-dH/(Rg*T))
    K = K0*np.exp(-dH_K/(Rg*T))
    q_star = q_s*b*partial_p**m/(1+b*partial_p**m) + K*partial_p
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

# #Fig 1a
# article_25 = np.load('25_graus.npy')
# article_40 = np.load('40_graus.npy') 
# article_60 = np.load('60_graus.npy') 

# plt.figure(figsize=(5, 5),dpi=200)
# T_interval = [25,40,60]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15
#     q_star_array = np.zeros((3,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,partial_p_interval[i])

#     plt.plot(partial_p_interval,q_star_array[0],label=f'{round(T-273.15)}°C',color=tol[j])

# plt.plot(article_25[:,0],article_25[:,1],color='black',linestyle=':',ms=3,marker='s',label='artigo')
# plt.plot(article_40[:,0],article_40[:,1],color='black',linestyle=':',ms=3,marker='s')
# plt.plot(article_60[:,0],article_60[:,1],color='black',linestyle=':',ms=3,marker='s')


# plt.xlabel("Partial pressure of CO$_2$, p (bar)")
# plt.ylabel(f"Loading (mmol/g {calf_type})")
# plt.xlim(0)
# plt.ylim(0,4)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_Fig1a.png',dpi=200)

# #Fig 1b
# plt.figure(figsize=(5, 5),dpi=200)
# T_interval = [25,40,60]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15
#     q_star_array = np.zeros((3,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,partial_p_interval[i])

#     plt.plot(partial_p_interval,q_star_array[1],label=f'{round(T-273.15)}°C',color=tol[j])

# plt.xlabel("Partial pressure of N$_2$, p (bar)")
# plt.ylabel(f"Loading (mmol/g {calf_type})")
# plt.xlim(0)
# plt.ylim(0,0.4)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_Fig1b.png',dpi=200)

# Fig 2a
# plt.figure(figsize=(5, 5),dpi=200)
# T_interval = [25,30,40,65,80]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15
#     q_star_array = np.zeros((3,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,partial_p_interval[i])

#     plt.plot(partial_p_interval,q_star_array[2],label=f'{round(T-273.15)}°C',color=tol[j])

# plt.xlabel("Partial pressure of H$_2$O, p (bar)")
# plt.ylabel(f"Loading (mmol/g {calf_type})")
# plt.xlim(0,0.45)
# plt.ylim(0,10)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_Fig2a.png',dpi=200)

#Fig 2 b
# plt.figure(figsize=(5, 5),dpi=200)
# T_interval = [120,150]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15
#     q_star_array = np.zeros((3,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,partial_p_interval[i])

#     plt.plot(partial_p_interval,q_star_array[2],label=f'{round(T-273.15)}°C',color=tol[j])

# plt.xlabel("Partial pressure of H$_2$O, p (bar)")
# plt.ylabel(f"Loading (mmol/g {calf_type})")
# plt.xlim(0,1)
# plt.ylim(0,8)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_Fig2b.png',dpi=200)


#Fig 4a

article_25 = np.load('25_graus.npy')
article_40 = np.load('40_graus.npy') 
article_60 = np.load('60_graus.npy') 
article_120 = np.load('120_graus.npy') 
article_150 = np.load('150_graus.npy') 

plt.figure(figsize=(5, 5),dpi=200)
T_interval = [25,40,60,120,150]
for j in range(len(T_interval)):
    T = T_interval[j] + 273.15
    q_star_array = np.zeros((3,N))
    for i in range(0,N):
        q_star_array[:,i] = calc_loading_single(T,partial_p_interval[i])

    plt.plot(partial_p_interval,q_star_array[0],label=f'{round(T-273.15)}°C',color=tol[j])

plt.plot(article_25[:,0],article_25[:,1],color='black',linestyle=':',ms=3,marker='s',label='artigo')
plt.plot(article_40[:,0],article_40[:,1],color='black',linestyle=':',ms=3,marker='s')
plt.plot(article_60[:,0],article_60[:,1],color='black',linestyle=':',ms=3,marker='s')
plt.plot(article_120[:,0],article_120[:,1],color='black',linestyle=':',ms=3,marker='s')
plt.plot(article_150[:,0],article_150[:,1],color='black',linestyle=':',ms=3,marker='s')

plt.xlabel("Partial pressure of CO$_2$, p (bar)")
plt.ylabel(f"Loading (mmol/g {calf_type})")
plt.xlim(0,1)
plt.ylim(0,4)
plt.legend()
plt.tight_layout()
plt.savefig('isoterma_Fig4a.png',dpi=200)

#Fig 4b
# plt.figure(figsize=(5, 5),dpi=200)
# T_interval = [25,40,60]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15
#     q_star_array = np.zeros((3,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,partial_p_interval[i])

#     plt.plot(partial_p_interval,q_star_array[1],label=f'{round(T-273.15)}°C',color=tol[j])

# plt.xlabel("Partial pressure of N$_2$, p (bar)")
# plt.ylabel(f"Loading (mmol/g {calf_type})")
# plt.xlim(0,1)
# plt.ylim(0,0.35)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_Fig4b.png',dpi=200)