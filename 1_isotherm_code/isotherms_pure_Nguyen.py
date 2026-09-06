import numpy as np
import matplotlib.pyplot as plt

R = 8.314                                   #J/mol.K
N = 100


#Parameters
#SSL
q_sb_ssl = np.array([3.927,5.658])          #mmol/g
b0_ssl = np.array([1.889e-8,8.143e-7])      #m³/mol   
dU_b_ssl = -np.array([41.19,17.96])*1e3     #kJ/mol         

#DSL
q_sb_dsl = np.array([2.387,2.387])          #mmol/g
b0_dsl = np.array([5.519e-7,8.139e-7])      #m³/mol   
dU_b_dsl = -np.array([35.06,17.96])*1e3     #kJ/mol         

q_sd = np.array([3.2711,3.271])         #mmol/g
d0 = np.array([5.187e-8,8.139e-7])      #m³/mol   
dU_d = -np.array([28.95,17.96])*1e3     #kJ/mol   


#Pure Component Isotherms
def calc_loading_single(T,conc,model):
    if model == 'SSL':
        q_sb = q_sb_ssl
        b0 = b0_ssl
        dU_b = dU_b_ssl
        b = b0*np.exp(-dU_b/(R*T))
        q_star = q_sb*b*conc/(1+b*conc)     
    elif model == 'DSL':
        q_sb = q_sb_dsl
        b0 = b0_dsl
        dU_b = dU_b_dsl
        b = b0*np.exp(-dU_b/(R*T))
        d = d0*np.exp(-dU_d/(R*T))
        q_star = q_sb*b*conc/(1+b*conc) + q_sd*d*conc/(1+d*conc)
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

# #Fig 2.4 Nguyen
# partial_p_interval = np.linspace(0,1,N)
# plt.figure(figsize=(5, 5),dpi=200)

# #SSL
# T_interval = [30,50]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15

#     concentration = partial_p_interval*1e5/(R*T)

#     q_star_array = np.zeros((2,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,concentration[i],model = 'SSL')

#     plt.plot(partial_p_interval,q_star_array[0],label=f'{T_interval[j]}°C',color=tol[j])

# #DSL
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15

#     concentration = partial_p_interval*1e5/(R*T)

#     q_star_array = np.zeros((2,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,concentration[i],model = 'DSL')

#     plt.plot(partial_p_interval,q_star_array[0],label=f'{T_interval[j]}°C',color=tol[j],ls='--')


# plt.xlabel("Partial pressure of CO$_2$, p (bar)")
# plt.ylabel(f"Loading CO$_2$ (mmol/g)")
# plt.xlim(0,1)
# plt.ylim(0,4)
# # plt.locator_params(axis='y', nbins=5)
# # plt.locator_params(axis='x', nbins=6)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_tese_Nguyen_Fig2.4.png',dpi=200)


#Fig 2.3b Nguyen
partial_p_interval = np.linspace(0,5,N)
plt.figure(figsize=(5, 5),dpi=200)

#SSL
T_interval = [30,40,50,60,80,100]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15

#     concentration = partial_p_interval*1e5/(R*T)

#     q_star_array = np.zeros((2,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,concentration[i],model = 'SSL')

#     plt.plot(partial_p_interval,q_star_array[0],label=f'{T_interval[j]}°C',color=tol[j])

# #DSL
for j in range(len(T_interval)):
    T = T_interval[j] + 273.15

    concentration = partial_p_interval*1e5/(R*T)

    q_star_array = np.zeros((2,N))
    for i in range(0,N):
        q_star_array[:,i] = calc_loading_single(T,concentration[i],model = 'DSL')

    plt.plot(partial_p_interval,q_star_array[0],label=f'{T_interval[j]}°C',color=tol[j],ls='-')


plt.xlabel("Partial pressure of CO$_2$, p (bar)")
plt.ylabel(f"Loading CO$_2$ (mmol/g)")
plt.xlim(0,5)
plt.ylim(0,4)
# plt.locator_params(axis='y', nbins=5)
# plt.locator_params(axis='x', nbins=6)
plt.legend()
plt.tight_layout()
plt.savefig('isoterma_tese_Nguyen_Fig2.3b.png',dpi=200)




# #Fig 2.3a Nguyen
# partial_p_interval = np.linspace(0,5,N)
# plt.figure(figsize=(5, 5),dpi=200)

# #SSL
# T_interval = [30,40,50,60,80,100]
# for j in range(len(T_interval)):
#     T = T_interval[j] + 273.15

#     concentration = partial_p_interval*1e5/(R*T)

#     q_star_array = np.zeros((2,N))
#     for i in range(0,N):
#         q_star_array[:,i] = calc_loading_single(T,concentration[i],model = 'SSL')

#     plt.plot(partial_p_interval,q_star_array[1],label=f'{T_interval[j]}°C',color=tol[j])

# plt.xlabel("Partial pressure of N$_2$, p (bar)")
# plt.ylabel(f"Loading N$_2$ (mmol/g)")
# plt.xlim(0,5)
# plt.ylim(0,1)
# # plt.locator_params(axis='y', nbins=5)
# # plt.locator_params(axis='x', nbins=6)
# plt.legend()
# plt.tight_layout()
# plt.savefig('isoterma_tese_Nguyen_Fig2.3a.png',dpi=200)