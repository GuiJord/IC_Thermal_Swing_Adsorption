import numpy as np
from scipy.integrate import solve_ivp

# Number of points
N = 30

save_solution_flag1 = True
print_iteration = True

# Cycle Phases
t_press = 15
t_ads = 30
t_blow = 30
t_evac = 40

lambda_press = 0.5
lambda_blow = 0.5
lambda_evac = 0.5

# Operating conditions for case studies
v_feed = 1 #m/s
T_feed = 298.15 #K
T_amb = 298.15
P_H = 1e5 #Pa

y_feed_CO2 = 0.15

# Constants=============================================================================================

R = 8.314

# Column Properties
L = 1 #m
r_in = 0.1445 #m
r_out = 0.162 #m
epsilon = 0.37
epsilon_prime = (1-epsilon)/epsilon
epsilon_p = 0.35
r_p = 1e-3 #m
d_p = 2*r_p
tortuosity = 3 
area = np.pi*r_in**2

# Properties and Constants
P_f = 1e5 #Pa
P_atm = 1e5
rho_s = 1130 #kg/m³
rho_w = 7800 #kg/m³
Cp_g = 30.7 #J/mol.K
Cp_ads = 30.7 #J/mol.K
Cp_s = 1070 #J/kg.K
Cp_w = 502 #J/kg.K

mu = 1.72*1e-5 #kg/m.s
D_m = 1.6*1e-5 #m²/s
gamma = 1.4
K_z = 0.09 #J/m.K.s
K_w = 16 #J/m.K.s
h_in = 0#8.6 #J/(m².K.s)
h_out = 0#2.5 #J/(m².K.s)

efficiency = 0.72

D_p = D_m/tortuosity


# Adsorption
b0 = np.array([8.65*1e-7,2.5*1e-6])
d0 = np.array([2.63*1e-8,0])

dU_b = np.array([-36641.21,-1.58*1e4])
dU_d = np.array([-35690.66,0])

q_sb = np.array([3.09,5.84])
q_sd = np.array([2.54,0])

# Constants for the dimensionless parameters
P0 = P_H
T0 = T_feed
v0 = v_feed
q_s0 = q_sb[1]

T_feed_ad = T_feed/T0
T_amb_ad = T_amb/T0

D_L = 0.7*D_m + 0.5*v0*d_p

Pe = v0*L/D_L

psi = R*T0*q_s0*rho_s/(P_H)*(1-epsilon)/epsilon

v_ad_feed = v_feed/v0


dz = L/N
dZ_ad = dz/L

dH_CO2 = dU_b[0] - R*T0
dH_N2 = dU_b[1] - R*T0

mole_0_cst = P0*v0/(R*T0)*epsilon*area
energy_0_cst = 1/efficiency*epsilon*area*v0*P0*gamma/(gamma-1)
energy_0_cst2 = (gamma-1)/gamma
P_f_ad = P_f/P0
P_atm_ad = P_atm/P0


def TVD_van_Leer(f,f_wall):
    f_0 = f[0] - 2*(f[0]-f_wall[0])
    r = (f[0] - f_0 + 1e-10)/(f[1]-f[0]+1e-10)
    phi = (r + np.abs(r))/(1 + np.abs(r))
    f_wall[1] = f[0] + 0.5*phi*(f[1]-f[0])

    for j in range(1,N-1):
        r = (f[j]-f[j-1]+1e-10)/(f[j+1]-f[j]+1e-10)
        phi = (r + np.abs(r))/(1 + np.abs(r))
        f_wall[j+1] = f[j] + 0.5*phi*(f[j+1]-f[j])

    return f_wall

def WENO(f,f_wall):
    f_0 = f[0] - 2*(f[0]-f_wall[0])    
    a0 = 2/3/(f[1]-f[0]+1e-10)**4
    a1 = 1/3/(f[0]-f_0+1e-10)**4
    den = a0 + a1
    f_wall[1] = a0/den*0.5*(f[0]+f[1]) + a1/den*0.5*(3*f[0]-f_0) 

    for j in range(1,N-1):
        a0 = 2/3/(f[j+1]-f[j]+1e-10)**4
        a1 = 1/3/(f[j]-f[j-1]+1e-10)**4
        den = a0 + a1
        f_wall[j+1] = a0/den*0.5*(f[j]+f[j+1]) + a1/den*0.5*(3*f[j]-f[j-1]) 

    return f_wall

def q_star_calc(concentration, temperature):
    b = b0*np.exp(-dU_b/(R*temperature))
    d = d0*np.exp(-dU_d/(R*temperature))

    A1 = q_sb*b/(1+np.sum(b*concentration))
    A2 = q_sd*d/(1+np.sum(d*concentration))

    q_star = concentration*(A1 + A2)
    return q_star


