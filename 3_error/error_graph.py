import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter, NullFormatter, FixedLocator
import scipy

#================================== Parameters - Start =======================================

tmax = 8000                     #s  Breakthrough maximum time
finite_volume_method = 'WENO'   #   Finite Volume Method ('WENO' or 'VANLEER')
solve_ivp_method = 'LSODA'      #   solve_ivp ODE solver ('LSODA' or 'BDF', etc)
adiabatic_operation = True      #   Adiabatic Column (TRUE or FALSE)

#Save progress ==============================================
save_progress_flag = False       #Save Progress (TRUE OR FALSE)
save_interval = 100
next_save = 0

#Operating Conditions ==============================================
P_feed = 1e5        #Pa     Pressure of feed    
T_feed = 298.15     #K      Temperature of feed
v_feed = 1          #m/s    Velocity of feed
T_amb = 298.15      #K      Ambient Temperature
P_H = 1e5           #Pa     High Pressure
y_feed_CO2 = 0.15   #       CO2 composition in the feed

#Constants ==============================================
R = 8.314           #J/mol.K    Ideal Gas Constant
P_atm = 1e5         #Pa         Atmospheric Pressure

#Column properties ==============================================
L = 1                   #m          Length
r_in = 0.1445           #m          Inside Radius
r_out = 0.162           #m          Outside Radius
epsilon = 0.37          #           Column Porosity
epsilon_prime = (1-epsilon)/epsilon
area = np.pi*r_in**2    #m²         Area of cross-section
rho_w = 7800            #kg/m³      Wall Density
Cp_w = 502              #J/kg.K     Wall Specific Heat
K_w = 16                #J/m.K.s    Wall Thermal Conductivity
if adiabatic_operation:
    h_in = 0            #J/(m².K.s) Heat Transfer Coefficient (inside) 
    h_out = 0           #J/(m².K.s) Heat Transfer Coefficient (outside)
else:
    h_in = 8.6          #J/(m².K.s) Heat Transfer Coefficient (inside) 
    h_out = 2.5         #J/(m².K.s) Heat Transfer Coefficient (outside)

#Particle properties ==============================================
epsilon_p = 0.35        #       Particle Porosity
r_p = 1e-3              #m      Radius      
d_p = 2*r_p             #m      Diameter
tortuosity = 3          #       Tortuosity
rho_s = 1130            #kg/m³  Solid Density
Cp_s = 1070             #J/kg.K Solid Specific Heat

#Fluid properties ==============================================
Cp_g = 30.7             #J/mol.K    Gas Specific Heat
mu = 1.72*1e-5          #kg/m.s     Fluid Viscosity
D_m = 1.6*1e-5          #m²/s       Molecular Diffusivity
D_p = D_m/tortuosity    #m²/s       Macropore Diffusivity
K_z = 0.09              #J/m.K.s    Gas Thermal Conductivity

#Adsorbed Phase properties ==============================================
Cp_ads = 30.7                           #J/mol.K    Adsorbed Phase Specific Heat
b0 = np.array([8.65*1e-7,2.5*1e-6])     #m³/mol     Isotherm Parameter
d0 = np.array([2.63*1e-8,0])            #m³/mol     Isotherm Parameter
dU_b = np.array([-36641.21,-1.58*1e4])  #J/mol      Activation Energy
dU_d = np.array([-35690.66,0])          #J/mol      Activation Energy
q_sb = np.array([3.09,5.84])            #mmol/g     Isotherm Parameter
q_sd = np.array([2.54,0])               #mmol/g     Isotherm Parameter
   
#Constants for the Dimensionless parameters ==============================================
P0 = P_H        #Pa         Pressure Dimensionless constant                            
T0 = T_feed     #K          Temperature Dimensionless constant
v0 = v_feed     #m/s        Velocity Dimensionless constant
q_s0 = q_sb[1]  #mmol/g     Loading Dimensionless constant

#Calculated Constants ==============================================
dH_CO2 = dU_b[0] - R*T0     #J/mol  Enthalpy of Adsorption of CO2
dH_N2 = dU_b[1] - R*T0      #J/mol  Enthalpy of Adsorption of N2
D_L = 0.7*D_m + 0.5*v0*d_p  #m²/s   Axial Dispersion
Pe = v0*L/D_L               #       Peclet Number
permeability = 4/150*r_p**2*(epsilon/(1-epsilon))**2    #m² Permeability (Darcy)
psi = R*T0*q_s0*rho_s/P0*(1-epsilon)/epsilon          #   constant (material balance)
den = rho_w*Cp_w*v0                                     #   constant (wall energy balance)
coeff1 = K_w/den/L                                      #   constant (wall energy balance)
coeff2 = 2*r_in*h_in/(r_out*r_out-r_in*r_in)*L/den      #   constant (wall energy balance)
coeff3 = 2*r_out*h_out/(r_out*r_out-r_in*r_in)*L/den    #   constant (wall energy balance)

#Dimensionless parameters ==============================================
P_feed_ad = P_feed/P0               #   Dimensionless Feed Pressure
T_feed_ad = T_feed/T0               #   Dimensionless Feed Temperature
v_ad_feed = v_feed/v0               #   Dimensionless Feed Velocity
P_atm_ad = P_atm/P0                 #   Dimensionless Ambient Pressure
T_amb_ad = T_amb/T0                 #   Dimensionless Ambient Temperature 
mu_ad = mu*v0/(P0*L)                #   Dimensionless viscosity
permeability_ad = permeability/L**2 #   Dimensionless Permeability

#Operation Constants ==============================================
efficiency = 0.72       #   Compression/evacuation efficiency
gamma = 1.4             #   Adiabatic constant

#Operation Dimensionless parameters ==============================================
mole_0_cst = P0*v0/(R*T0)*epsilon*area
energy_0_cst = 1/efficiency*epsilon*area*v0*P0*gamma/(gamma-1)
energy_0_cst2 = (gamma-1)/gamma
#======================================== Constants - End

#================================== Parameters - End =========================================

N_list = np.array([10,20,30,40,50,60,70,80,90,100,200,300,400,500])
# error_MB = np.array([1.403560e-03,1.314498e-03,7.666741e-04,6.799225e-04,7.002546e-04,1.176237e-03,7.661735e-04,9.696675e-04,7.689534e-04,9.555427e-04,9.376370e-04,8.592614e-04,9.809059e-04,9.602344e-04])
# error_EB = np.array([6.966357e-03,3.428899e-03,2.310138e-03,1.721447e-03,1.363240e-03,1.139398e-03,9.424263e-04,7.998989e-04,7.247072e-04,6.247355e-04,2.796099e-04,1.805532e-04,1.265069e-04,9.255491e-05])
# time_list = np.array([2.459e+01,1.052e+02,2.666e+02,4.443e+02,8.448e+02,1.727e+03,1.973e+03,2.591e+03,3.498e+03,4.529e+03,1.146e+04,2.049e+04,2.630e+04,3.754e+04])

error_MB_list = np.array([])
error_EB_list = np.array([])


for N in N_list:

    dz = L/N                    #m      Z-axis discretization dz
    dZ_ad = dz/L                        #   Dimensionless dZ 
    dPdz_aux = -permeability_ad/mu_ad/dZ_ad    #auxiliary variable    v = aux*dP

    t = np.load(f'./run_{N}/t.npy')
    P_ad = np.load(f'./run_{N}/P_ad.npy')
    T_ad = np.load(f'./run_{N}/T_ad.npy')
    v_ad = np.load(f'./run_{N}/v_ad.npy')
    y_CO2 = np.load(f'./run_{N}/y_CO2.npy')
    x_CO2 = np.load(f'./run_{N}/x_CO2.npy')
    x_N2 = np.load(f'./run_{N}/x_N2.npy')
    Tw_ad = np.load(f'./run_{N}/Tw_ad.npy')

    #================================== Error Mesure - Start =======================================

    Z = np.linspace(0,1,N)
    #Mass Balance Error ==============================================
    #Breakthrough Boundary Conditions ==============================================
    P_ad_in        = P_ad[0,:] - v_ad_feed/dPdz_aux/2
    P_ad_out       = 1
    v_ad_in        = v_ad_feed
    v_ad_out       = 2*dPdz_aux*(P_ad_out-P_ad[-1,:])
    y_ad_CO2_in    = (y_CO2[0,:]+y_feed_CO2*v_ad_in*Pe*dZ_ad*0.5)/(1+v_ad_in*Pe*dZ_ad*0.5)
    y_ad_CO2_out   = y_CO2[-1,:]
    rho_g          = P0*P_ad[0,:]/(R*T0*T_ad[0,:])
    Pe_h           = epsilon*v0*L/K_z*rho_g*Cp_g
    T_ad_in        = (T_ad[0,:]+v_ad_in*Pe_h*dZ_ad*0.5)/(1+v_ad_in*Pe_h*dZ_ad*0.5)
    T_ad_out       = T_ad[-1,:]

    # P_ad_in        = P_ad[0,:]
    # P_ad_out       = P_ad[-1,:]
    # v_ad_in        = v_ad_feed
    # v_ad_out       = v_ad[:]
    # y_ad_CO2_in    = y_CO2[0,:]
    # y_ad_CO2_out   = y_CO2[-1,:]
    # rho_g          = P0*P_ad[0,:]/(R*T0*T_ad[0,:])
    # Pe_h           = epsilon*v0*L/K_z*rho_g*Cp_g
    # T_ad_in        = T_ad[0,:]
    # T_ad_out       = T_ad[-1,:]

    #Total Mass ==============================================
    integrand = v_ad_in*P_ad_in/T_ad_in
    n_in = P0/(R*T0)*epsilon*area*L*scipy.integrate.simpson(integrand,x=t)
    integrand = v_ad_out*P_ad_out/T_ad_out
    n_out = P0/(R*T0)*epsilon*area*L*scipy.integrate.simpson(integrand,x=t)
    #CO2 ==============================================
    #Inlet
    integrand = y_ad_CO2_in*v_ad_in*P_ad_in/T_ad_in
    n_CO2_in = P0/(R*T0)*epsilon*area*L*scipy.integrate.simpson(integrand,x=t)
    #Outlet
    integrand = y_ad_CO2_out*v_ad_out*P_ad_out/T_ad_out
    n_CO2_out = P0/(R*T0)*epsilon*area*L*scipy.integrate.simpson(integrand,x=t)
    #Accumulated in gas phase
    integrand = y_CO2[:,-1]*P_ad[:,-1]/T_ad[:,-1] - y_CO2[:,0]*P_ad[:,0]/T_ad[:,0]
    n_CO2_acc_g = P0/(R*T0)*epsilon*area*L*scipy.integrate.simpson(integrand,x=Z)
    #Accumulated in adsorbed phase
    integrand = x_CO2[:,-1] - x_CO2[:,0]
    n_CO2_acc_ads = q_s0*rho_s*(1-epsilon)*area*L*scipy.integrate.simpson(integrand,x=Z)
    #Accumulated Total
    n_CO2_acc = n_CO2_acc_g + n_CO2_acc_ads
    #Accumulated N2 in adsorbed phase
    integrand = x_N2[:,-1] - x_N2[:,0]
    n_N2_acc_ads = q_s0*rho_s*(1-epsilon)*area*L*scipy.integrate.simpson(integrand,x=Z)
    #Accumulated total in gas phase
    integrand = P_ad[:,-1]/T_ad[:,-1]
    n_acc_g = P0/(R*T0)*epsilon*area*L*scipy.integrate.simpson(integrand,x=Z)

    #Error Mass Balance ==============================================
    # error_MB = np.abs(n_CO2_in-n_CO2_out-n_CO2_acc)/np.abs(n_CO2_in)
    error_MB = np.abs(n_CO2_in-n_CO2_out-n_CO2_acc)/n_CO2_acc

    #======================================== Mass Balance Error - End

    #Energy Balance Error ==============================================
    rho_g = P_ad*P0/R/T_ad/T0
    #Heat Inlet ==============================================
    integrand = v_ad_in*T_ad_in*rho_g[0,:]
    heat_in = epsilon*area*L*T0*Cp_g*scipy.integrate.simpson(integrand,x=t)
    #Heat Outlet ==============================================
    integrand = v_ad_out*T_ad_out*rho_g[-1,:]
    heat_out = epsilon*area*L*T0*Cp_g*scipy.integrate.simpson(integrand,x=t)
    #Heat Generated ==============================================
    integrand = -(dH_CO2*(x_CO2[:,-1]-x_CO2[:,0]) + dH_N2*(x_N2[:,-1]-x_N2[:,0]))
    heat_gen = (1-epsilon)*area*L*q_s0*rho_s*scipy.integrate.simpson(integrand,x=Z)
    #Heat Accumulated in solid ==============================================
    integrand = T_ad[:,-1] - T_ad[:,0]
    heat_acc_s = (1-epsilon)*area*L*T0*Cp_s*rho_s*scipy.integrate.simpson(integrand,x=Z)
    #Heat Accumulated in Fluid ==============================================
    integrand = T_ad[:,-1]*rho_g[:,-1] - T_ad[:,0]*rho_g[:,0]
    heat_acc_f = epsilon*area*L*T0*Cp_g*scipy.integrate.simpson(integrand,x=Z)
    #Heat Accumulated in Adsorbed ==============================================
    integrand = (T_ad[:,-1]*x_CO2[:,-1] - T_ad[:,0]*x_CO2[:,0])+(T_ad[:,-1]*x_N2[:,-1] - T_ad[:,0]*x_N2[:,0])
    heat_acc_ads = (1-epsilon)*area*L*T0*Cp_ads*q_s0*rho_s*scipy.integrate.simpson(integrand,x=Z)
    #Error Energy Balance ==============================================
    # error_EB = np.abs(heat_in + heat_gen - (heat_out + heat_acc_s + heat_acc_f + heat_acc_ads))/np.abs(heat_in)
    error_EB = np.abs(heat_in + heat_gen - (heat_out + heat_acc_s + heat_acc_f + heat_acc_ads))/(heat_out + heat_acc_s + heat_acc_f + heat_acc_ads)
    #======================================== Energy Balance Error - End

    error_MB_list = np.append(error_MB_list,error_MB)
    error_EB_list = np.append(error_EB_list,error_EB)
    #================================== Error Mesure - End =========================================


MB_arvind = np.load('MB.npy')
EB_arvind = np.load('EB.npy')


def error_graph():
    fig, ax1 = plt.subplots(figsize=(9,6), dpi=200)
    ax2 = ax1.twinx()

    # Error MB
    ax1.plot(
        N_list,
        error_MB_list*100,
        color='green',
        linestyle='-',
        ms=6,marker='s',
        label='Erro BM'
    )

    ax1.plot(
        MB_arvind[:,0],
        MB_arvind[:,1],
        color='black',
        linestyle='-',
        ms=6,marker='s',
        label='Erro BM artigo'
    )

    # Error EB
    ax2.plot(
        N_list,
        error_EB_list*100,
        color='red',
        linestyle='-',
        ms=6,marker='^',
        label='Erro BE'
    )

    ax2.plot(
        EB_arvind[:,0],
        EB_arvind[:,1],
        color='black',
        linestyle='-',
        ms=6,marker='^',
        label='Erro BE artigo'
    )

    ax1.set_box_aspect(6/9)
    # ax1.set_ylim(0.01,6)
    # ax2.set_ylim(0.006,3)
    ax1.set_xlim(9,600)

    # Labels
    ax1.set_xlabel("Número de Volumes Finitos [-]", fontsize=13)
    ax1.set_ylabel("Erro Balanço de Massa [%]", fontsize=13)
    ax2.set_ylabel("Erro Balanço de Energia [%]", fontsize=13)

    # Log scales
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')

    # Tick appearance
    ax1.tick_params(axis='both', which='major',
                    direction='in', length=7, width=1.2,
                    labelsize=12, top=True)

    ax1.tick_params(axis='both', which='minor',
                    direction='in', length=4, width=0.8,
                    top=True)

    ax2.tick_params(axis='y', which='major',
                    direction='in', length=7, width=1.2,
                    labelsize=12, right=True)

    ax2.tick_params(axis='y', which='minor',
                    direction='in', length=4, width=0.8,
                    right=True)

    for x in [10, 20, 30, 50, 75, 100, 200, 300]:
        ax1.axvline(x,
                    color='gray',
                    linestyle='--',
                    linewidth=0.8,
                    alpha=0.3,
                    zorder=0)

    # Grid
    ax1.grid(which='major', linestyle='--', alpha=0.3)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=11)

    plt.tight_layout()
    plt.savefig("error_graph.png", dpi=200)

error_graph()