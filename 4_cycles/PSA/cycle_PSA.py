import numpy as np
import scipy
import time
import matplotlib.pyplot as plt

#================================== Parameters - Start =======================================

N = 30                          #   Number of Volumes
tmax = 8000                     #s  Breakthrough maximum time
finite_volume_method = 'WENO'   #   Finite Volume Method ('WENO' or 'VANLEER')
solve_ivp_method = 'LSODA'      #   solve_ivp ODE solver ('LSODA' or 'BDF', etc)
adiabatic_operation = False      #   Adiabatic Column (TRUE or FALSE)
operation_step = 'none'         #   Operation Step

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
P_L = 0.9e5
P_L_final = 0.1e5
P_I = 0.95e5
P_I_final = 0.2e5
lambda_press = 0.5
lambda_blow = 0.5
lambda_evac = 0.5
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
P0 = 1e5        #Pa         Pressure Dimensionless constant                            
T0 = T_feed     #K          Temperature Dimensionless constant
v0 = v_feed     #m/s        Velocity Dimensionless constant
q_s0 = q_sb[1]  #mmol/g     Loading Dimensionless constant

#Calculated Constants ==============================================
dz = L/N                    #m      Z-axis discretization dz
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
dZ_ad = dz/L                        #   Dimensionless dZ 
mu_ad = mu*v0/(P0*L)                #   Dimensionless viscosity
permeability_ad = permeability/L**2 #   Dimensionless Permeability
dPdz_aux = -permeability_ad/mu_ad/dZ_ad    #auxiliary variable    v = aux*dP

#Operation Constants ==============================================
efficiency = 0.72       #   Compression/evacuation efficiency
gamma = 1.4             #   Adiabatic constant

#Operation Dimensionless parameters ==============================================
mole_0_cst = P0*v0/(R*T0)*epsilon*area
energy_0_cst = 1/efficiency*epsilon*area*v0*P0*gamma/(gamma-1)
energy_0_cst2 = (gamma-1)/gamma
#======================================== Constants - End

#================================== Parameters - End =========================================


#================================== Numerical Methods (Finite Volumes) - Start =======================================

#TVD van Leer ==============================================
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

#WENO ==============================================
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

#Method Selection ==============================================
if finite_volume_method == 'WENO':
    finite_volume = WENO
elif finite_volume_method == 'VANLEER':
    finite_volume = TVD_van_Leer
else:
    raise ValueError("Choose Finite Volume Method")
#================================== Numerical Methods (Finite Volumes) - End =========================================


#================================== Isotherm Calculation - Start =======================================

#Loading ==============================================
def q_star_calc(concentration, temperature):
    b = b0*np.exp(-dU_b/(R*temperature))
    d = d0*np.exp(-dU_d/(R*temperature))

    A1 = q_sb*b/(1+np.sum(b*concentration))
    A2 = q_sd*d/(1+np.sum(d*concentration))

    q_star = concentration*(A1 + A2)
    return q_star

#Solid phase LDF constant and equilibrium loading ==============================================
def solid_phase_calc(concentration,temperature):
    b = b0*np.exp(-dU_b/(R*temperature))
    d = d0*np.exp(-dU_d/(R*temperature))

    A1 = q_sb*b/(1+np.sum(b*concentration))
    A2 = q_sd*d/(1+np.sum(d*concentration))

    q_star_over_concentration = A1 + A2
    alpha = 1/q_star_over_concentration/rho_s*(15*epsilon_p*D_p/(r_p*r_p))*L/v0
    x_star = q_star_over_concentration*concentration/q_s0
    return alpha, x_star

#================================== Isotherm Calculation - End =========================================


#================================== ODE System - Start =======================================

def ode_system(t, variable_vec):
    global operation_step,t_actual
    global next_save, t_global, P_ad_global, T_ad_global, y_CO2_global, x_CO2_global, x_N2_global, Tw_ad_global
     
    #Arrays ==============================================
    # ODE Variables
    P_ad,T_ad,y_CO2,x_CO2,x_N2,Tw_ad = np.split(variable_vec, 6)
    # Walls
    P_ad_wall = np.zeros(N+1)
    T_ad_wall = np.zeros(N+1)
    y_CO2_wall = np.zeros(N+1)
    v_ad_wall = np.zeros(N+1)
    # Gradients
    dTdz_forward = np.zeros(N)
    dTdz_backward = np.zeros(N)
    dydz_CO2_forward = np.zeros(N)
    dydz_CO2_backward = np.zeros(N)
    dTwdz_forward = np.zeros(N)
    dTwdz_backward = np.zeros(N)
    # Differential Equations
    dP_dt = np.zeros(N)
    dT_dt = np.zeros(N)
    dy_dt_CO2 = np.zeros(N)
    dx_dt_CO2 = np.zeros(N)
    dx_dt_N2 = np.zeros(N)
    dTw_dt = np.zeros(N)
    #======================================== Arrays - End
    
    #Boundary Conditions ==============================================
    if operation_step == 'pressurization':
        #Pressure
        P_ad_wall[0]    = 1/P0*(P_H-(P_H-P_L)*np.exp(-lambda_press*t*L/v0))
        P_ad_wall[-1]   = P_ad[-1]
        #Velocity (Limiter for inlet)
        # v_ad_wall[0]    = v_ad_feed*max((P_H - P_ad_wall[0]*P0)/(P_H - P_L), 0.0)
        v_ad_wall[0]   = -2*dPdz_aux*(P_ad_wall[0]-P_ad[0])
        v_ad_wall[-1]   = 2*dPdz_aux*(P_ad_wall[-1]-P_ad[-1])
        #Composition
        y_CO2_wall[0]   = (y_CO2[0]+y_feed_CO2*v_ad_wall[0]*Pe*dZ_ad*0.5)/(1+v_ad_wall[0]*Pe*dZ_ad*0.5)
        y_CO2_wall[-1]  = y_CO2[-1]
        #Temperature
        rho_g           = P0*P_ad[0]/(R*T0*T_ad[0])
        Pe_h            = epsilon*v0*L/K_z*rho_g*Cp_g
        T_ad_wall[0]    = (T_ad[0]+v_ad_wall[0]*Pe_h*dZ_ad*0.5)/(1+v_ad_wall[0]*Pe_h*dZ_ad*0.5)
        T_ad_wall[-1]   = T_ad[-1]
        #Wall Temperature
        Tw_ad[0]        = T_amb_ad
        Tw_ad[-1]       = T_amb_ad

        #Spatial Derivatives
        dydz_CO2_backward[0]    = -v_ad_wall[0]*Pe*(y_feed_CO2-y_CO2_wall[0])
        dydz_CO2_forward[-1]    = 0
        
        dTdz_backward[0]        = -v_ad_wall[0]*Pe_h*(T_feed_ad-T_ad_wall[0])
        dTdz_forward[-1]        = 0

        dTwdz_backward[0]       = 0
        dTwdz_forward[-1]       = 0
    
    elif operation_step == 'adsorption':
        #Pressure
        P_ad_wall[0]    = P_ad[0] - v_ad_feed/dPdz_aux/2
        P_ad_wall[-1]   = P_H/P0
        #Velocity
        v_ad_wall[0]    = v_ad_feed
        v_ad_wall[-1]   = 2*dPdz_aux*(P_ad_wall[-1]-P_ad[-1])
        #Composition
        y_CO2_wall[0]   = (y_CO2[0]+y_feed_CO2*v_ad_feed*Pe*dZ_ad*0.5)/(1+v_ad_feed*Pe*dZ_ad*0.5)
        y_CO2_wall[-1]  = y_CO2[-1]
        #Temperature
        rho_g           = P0*P_ad[0]/(R*T0*T_ad[0])
        Pe_h            = epsilon*v0*L/K_z*rho_g*Cp_g
        T_ad_wall[0]    = (T_ad[0]+v_ad_feed*Pe_h*dZ_ad*0.5)/(1+v_ad_feed*Pe_h*dZ_ad*0.5)
        T_ad_wall[-1]   = T_ad[-1]
        #Wall Temperature
        Tw_ad[0]        = T_amb_ad
        Tw_ad[-1]       = T_amb_ad

        #Spatial Derivatives
        dydz_CO2_backward[0]    = -v_ad_wall[0]*Pe*(y_feed_CO2-y_CO2_wall[0])
        dydz_CO2_forward[-1]    = 0
        
        dTdz_backward[0]        = -v_ad_wall[0]*Pe_h*(T_feed_ad-T_ad_wall[0])
        dTdz_forward[-1]        = 0

        dTwdz_backward[0]       = 0
        dTwdz_forward[-1]       = 0
    elif operation_step == 'blowdown':
        #Pressure
        P_ad_wall[0] = P_ad[0]
        P_ad_wall[-1] = 1/P0*(P_I+(P_H-P_I)*np.exp(-lambda_blow*t*L/v0))
        #Velocity
        v_ad_wall[0] = 0
        v_ad_wall[-1] = 2*dPdz_aux*(P_ad_wall[-1]-P_ad[-1])
        # v_ad_wall[-1] = v_ad_feed*max((P_H - P_ad_wall[-1]*P0)/(P_H - P_I), 0.0)
        # print('P_w:',P_ad_wall[-1],'P',P_ad[-1],'v:',v_ad_wall[-1])
        #Composition
        y_CO2_wall[0] = y_CO2[0]
        y_CO2_wall[-1] = y_CO2[-1]
        #Temperature
        T_ad_wall[0] = T_ad[0]
        T_ad_wall[-1] = T_ad[-1]
        #Wall Temperature
        Tw_ad[0] = T_amb_ad
        Tw_ad[-1] = T_amb_ad

        #Spatial Derivatives
        dydz_CO2_backward[0]   = 0
        dydz_CO2_forward[-1]   = 0
        
        dTdz_backward[0]       = 0
        dTdz_forward[-1]       = 0
    
        dTwdz_backward[0]  = 0
        dTwdz_forward[-1]  = 0

    elif operation_step == 'evacuation':
        P_ad_wall[0] = 1/P0*(P_L+(P_I-P_L)*np.exp(-lambda_evac*t*L/v0))
        P_ad_wall[-1] = P_ad[-1]
        
        v_ad_wall[0]   = 2*dPdz_aux*(P_ad[0]-P_ad_wall[0])
        # v_ad_wall[0] = -v_ad_feed*max((P_L - P_ad_wall[0]*P0)/(P_L - P_I), 0.0)
        v_ad_wall[-1] = 0
        # print('P_w:',P_ad_wall[0],'P',P_ad[0],'v:',v_ad_wall[0])

        y_CO2_wall[0] = y_CO2[0]
        y_CO2_wall[-1] = y_CO2[-1]

        T_ad_wall[0] = T_ad[0]
        T_ad_wall[-1] = T_ad[-1]

        Tw_ad[0] = T_amb_ad
        Tw_ad[-1] = T_amb_ad

        dydz_CO2_backward[0]   = 0
        dydz_CO2_forward[-1]   = 0

        dTdz_backward[0]       = 0
        dTdz_forward[-1]       = 0
        
        dTwdz_backward[0]  = 0
        dTwdz_forward[-1]  = 0
    
    #======================================== Boundary Conditions - End    
    
    #Finite Volumes ==============================================
    P_ad_wall = finite_volume(P_ad, P_ad_wall)
    T_ad_wall = finite_volume(T_ad, T_ad_wall)
    y_CO2_wall = finite_volume(y_CO2, y_CO2_wall)
    #======================================== Finite Volumes - End

    #Velocity Calculation ==============================================
    for j in range(0,N-1):
        v_ad_wall[j+1] = dPdz_aux*(P_ad[j+1]-P_ad[j])   #Darcy Equation
    #======================================== Velocity Calculation - End
    
    #Spatial Derivatives ==============================================
    #forward
    for j in range(N-1):
        dTdz_forward[j] = (T_ad[j+1]-T_ad[j])/dZ_ad
        dydz_CO2_forward[j] = (y_CO2[j+1]-y_CO2[j])/dZ_ad
        dTwdz_forward[j] = (Tw_ad[j+1]-Tw_ad[j])/dZ_ad
    #backward
    for j in range(1,N):
        dTdz_backward[j] = (T_ad[j]-T_ad[j-1])/dZ_ad
        dydz_CO2_backward[j] = (y_CO2[j]-y_CO2[j-1])/dZ_ad
        dTwdz_backward[j] = (Tw_ad[j]-Tw_ad[j-1])/dZ_ad
    #======================================== Spatial Derivatives - End
        
    #Balance Equations (Time Derivatives) ==============================================
    for j in range(N):
        #Auxiliary variables
        y_CO2_half_p = y_CO2_wall[j+1]
        P_half_p = P_ad_wall[j+1]
        T_half_p = T_ad_wall[j+1]
        v_half_p = v_ad_wall[j+1]
        y_CO2_half_n = y_CO2_wall[j]
        P_half_n = P_ad_wall[j]
        T_half_n = T_ad_wall[j]
        v_half_n = v_ad_wall[j]

        #Solid Phase Balance ==============================================
        concentration = np.array([y_CO2[j],1-y_CO2[j]])*(P_ad[j]*P0)/(R*T_ad[j]*T0)
        alpha, x_star = solid_phase_calc(concentration,T_ad[j]*T0)

        dx_dt_CO2[j] = alpha[0]*(x_star[0] - x_CO2[j])
        dx_dt_N2[j] = alpha[1]*(x_star[1] - x_N2[j])
        #Total Mass and Column Energy Balance ==============================================
        dP_dt_prime = -T_ad[j]/dZ_ad*(P_half_p/T_half_p*v_half_p - P_half_n/T_half_n*v_half_n)\
                        -psi*T_ad[j]*(dx_dt_CO2[j]+dx_dt_N2[j])
        den = rho_s*Cp_s + q_s0*rho_s*Cp_ads*(x_CO2[j]+x_N2[j])
        Omega1 = K_z/(v0*epsilon*L)/(epsilon_prime*den)
        Omega2 = Cp_g*P0/(R*T0)/(epsilon_prime*den)
        Omega3 = Cp_ads*q_s0*rho_s/den
        Omega4 = 2*h_in*L/(r_in*v0)/((1-epsilon)*den)
        sigma_CO2 = q_s0*rho_s/T0*(-dH_CO2)/den
        sigma_N2 = q_s0*rho_s/T0*(-dH_N2)/den

        dT_dt[j] = (Omega1/dZ_ad*(dTdz_forward[j]-dTdz_backward[j])\
                    -Omega2/dZ_ad*(v_half_p*P_half_p-v_half_n*P_half_n)\
                    -Omega3*T_ad[j]*(dx_dt_CO2[j]+dx_dt_N2[j])\
                    +(sigma_CO2*dx_dt_CO2[j]+sigma_N2*dx_dt_N2[j])\
                    -Omega4*(T_ad[j]-Tw_ad[j])\
                    -Omega2*dP_dt_prime)/(1+Omega2*P_ad[j]/T_ad[j])
        
        dP_dt[j] = dP_dt_prime + P_ad[j]/T_ad[j]*dT_dt[j]
        #Wall Energy Balance ==============================================
        dTw_dt[j] = coeff1/dZ_ad*(dTwdz_forward[j]-dTwdz_backward[j])\
                        +coeff2*(T_ad[j]-Tw_ad[j])\
                        -coeff3*(Tw_ad[j]-T_amb_ad)
        #Component Mass Balance ==============================================
        dy_dt_CO2[j] = 1/Pe*T_ad[j]/P_ad[j]/dZ_ad*(P_half_p/T_half_p*dydz_CO2_forward[j] - P_half_n/T_half_n*dydz_CO2_backward[j])\
                    -T_ad[j]/P_ad[j]/dZ_ad*(y_CO2_half_p*P_half_p/T_half_p*v_half_p - y_CO2_half_n*P_half_n/T_half_n*v_half_n)\
                    -psi*T_ad[j]/P_ad[j]*dx_dt_CO2[j]\
                    -y_CO2[j]/P_ad[j]*dP_dt[j]\
                    +y_CO2[j]/T_ad[j]*dT_dt[j]
    #======================================== Balance Equations (Time Derivatives) - End
    
    #Save progress ==============================================
    if save_progress_flag and t > next_save:
        t_global = np.append(t_global,t)
        P_ad_global = np.vstack((P_ad_global,P_ad))
        T_ad_global = np.vstack((T_ad_global,T_ad))
        y_CO2_global = np.vstack((y_CO2_global,y_CO2))
        x_CO2_global = np.vstack((x_CO2_global,x_CO2))
        x_N2_global = np.vstack((x_N2_global,x_N2))
        Tw_ad_global = np.vstack((Tw_ad_global,Tw_ad))

        np.save("t",t_global.T)
        np.save("P_ad",P_ad_global.T)
        np.save("T_ad",T_ad_global.T)
        np.save("y_CO2",y_CO2_global.T)
        np.save("x_CO2",x_CO2_global.T)
        np.save("x_N2",x_N2_global.T)
        np.save("Tw_ad",Tw_ad_global.T)
        next_save += save_interval 
    #======================================== Save progress - End
    
    #Derivative array ==============================================
    derivative_vec = np.concatenate([dP_dt,dT_dt,dy_dt_CO2,dx_dt_CO2,dx_dt_N2,dTw_dt])

    #Progress ==============================================
    print(f'Progress:\t{t/(t_step*v0/L)*100:3.1f}%',end="\r")
    return derivative_vec

#================================== ODE System - End =========================================


#================================== Result Mesure - Start =======================================

def result_calc(t,variable_vec):
    global operation_step
    P_ad,T_ad,y_CO2,x_CO2,x_N2,Tw_ad = np.split(variable_vec, 6)

    Z = np.linspace(0,1,N)
    #Mass Balance Error ==============================================
    if operation_step == 'pressurization':
        #Pressurization Boundary Conditions ==============================================
        P_ad_in         = 1/P0*(P_H-(P_H-P_L)*np.exp(-lambda_press*t*L/v0))#np.array([1/P0*(P_H-(P_H-P_L)*np.exp(-lambda_press*t_*L/v0)) for t_ in t])#P_ad[0,:]#
        P_ad_out        = P_ad[-1,:]
        v_ad_in         = -2*dPdz_aux*(P_ad_in-P_ad[0,:])#v_ad_feed*np.array([max((P_H - P_ad_in_*P0)/(P_H - P_L), 0.0) for P_ad_in_ in P_ad_in])#dPdz_aux*(P_ad[1,:]-P_ad[0,:])#
        v_ad_out        = 2*dPdz_aux*(P_ad_out-P_ad[-1,:])
        y_ad_CO2_in     = (y_CO2[0,:]+y_feed_CO2*v_ad_in*Pe*dZ_ad*0.5)/(1+v_ad_in*Pe*dZ_ad*0.5)
        y_ad_CO2_out    = y_CO2[-1,:]
        rho_g           = P0*P_ad[0,:]/(R*T0*T_ad[0,:])
        Pe_h            = epsilon*v0*L/K_z*rho_g*Cp_g
        T_ad_in         = (T_ad[0,:]+v_ad_in*Pe_h*dZ_ad*0.5)/(1+v_ad_in*Pe_h*dZ_ad*0.5)
        T_ad_out        = T_ad[-1,:]
    
    elif operation_step == 'adsorption':
        #Adsorption Boundary Conditions ==============================================
        P_ad_in         = P_ad[0,:] - v_ad_feed/dPdz_aux/2
        P_ad_out        = P_H/P0
        v_ad_in         = v_ad_feed
        v_ad_out        = 2*dPdz_aux*(P_ad_out-P_ad[-1,:])
        y_ad_CO2_in     = (y_CO2[0,:]+y_feed_CO2*v_ad_in*Pe*dZ_ad*0.5)/(1+v_ad_in*Pe*dZ_ad*0.5)
        y_ad_CO2_out    = y_CO2[-1,:]
        rho_g           = P0*P_ad[0,:]/(R*T0*T_ad[0,:])
        Pe_h            = epsilon*v0*L/K_z*rho_g*Cp_g
        T_ad_in         = (T_ad[0,:]+v_ad_in*Pe_h*dZ_ad*0.5)/(1+v_ad_in*Pe_h*dZ_ad*0.5)
        T_ad_out        = T_ad[-1,:]

    elif operation_step == 'blowdown':
        P_ad_in         = P_ad[0,:]
        P_ad_out        = 1/P0*(P_I+(P_H-P_I)*np.exp(-lambda_blow*t*L/v0))
        v_ad_in         = 0
        v_ad_out        = 2*dPdz_aux*(P_ad_out-P_ad[-1,:])
        y_ad_CO2_in     = y_CO2[0,:]
        y_ad_CO2_out    = y_CO2[-1,:]
        T_ad_in         = T_ad[0,:]
        T_ad_out        = T_ad[-1,:]

        plt.figure(figsize=(9, 4))
        plt.plot(t, v_ad_out,color='dodgerblue')
        plt.ylabel('v$_{out}$ [-]')
        plt.xlabel('Time [-]')
        plt.xlim(left=0)
        plt.tight_layout()
        plt.savefig('velocity_time.png',dpi=200)

    elif operation_step == 'evacuation':
        P_ad_out        = 1/P0*(P_L+(P_I-P_L)*np.exp(-lambda_evac*t*L/v0))
        P_ad_in         = P_ad[-1,:]
        v_ad_out        = np.abs(2*dPdz_aux*(P_ad[0,:]-P_ad_out))
        v_ad_in         = 0
        y_ad_CO2_out    = y_CO2[0,:]
        y_ad_CO2_in     = y_CO2[-1,:]
        T_ad_out        = T_ad[0,:]
        T_ad_in         = T_ad[-1,:]

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
    error_MB = np.abs(n_CO2_in-n_CO2_out-n_CO2_acc)/np.abs(n_CO2_acc)
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
    heat_lost = heat_out + heat_acc_s + heat_acc_f + heat_acc_ads
    error_EB = np.abs(heat_in + heat_gen - (heat_out + heat_acc_s + heat_acc_f + heat_acc_ads))/np.abs(heat_out + heat_acc_s + heat_acc_f + heat_acc_ads)
    #======================================== Energy Balance Error - End

    #Print result ==============================================
    print(f'\tMass Balance Error:\t{error_MB*100:.2f}%')
    print(f'\tEnergy Balance Error:\t{error_EB*100:.2f}%')
    #Result Important data ==============================================
    result_data = np.array([n_CO2_in,n_CO2_out,n_CO2_acc,n_in,n_out,heat_in,heat_gen,heat_lost])
    return result_data

#================================== Result Mesure - End =========================================





def adsorption(variable_vec):
    global operation_step
    operation_step = 'adsorption'

    time_interval = np.array([0,t_ads])*v0/L

    #Solve_IVP ==============================================
    start = time.time()
    sol = scipy.integrate.solve_ivp(ode_system, time_interval, variable_vec, method=solve_ivp_method, rtol=1e-8, atol=1e-10)
    end = time.time()
    print(f"\nProcessing Time: {end - start:.2f} s")

    #Results ==============================================
    t = sol.t
    variable_vec = sol.y
    P_ad,T_ad,y_CO2,x_CO2,x_N2,Tw_ad = np.split(variable_vec, 6)

    #Velocity result in outlet ==============================================
    N_time_steps = len(t)
    v_ad = np.zeros(N_time_steps)
    for i in range(N_time_steps):
        P_ad_wall = np.zeros(N+1)
        P_ad_wall = WENO(P_ad[:,i],P_ad_wall)
        P_ad_wall[0] = P_ad[0,i] - v_ad_feed/dPdz_aux/2
        P_ad_wall[-1] = 1
        v_ad[i] = dPdz_aux*(P_ad_wall[-1]-P_ad_wall[-2])

    #Save ==============================================
    np.save(f"t_{operation_step}",t)
    np.save(f"P_ad_{operation_step}",P_ad)
    np.save(f"T_ad_{operation_step}",T_ad)
    np.save(f"v_ad_{operation_step}",v_ad)
    np.save(f"y_CO2_{operation_step}",y_CO2)
    np.save(f"x_CO2_{operation_step}",x_CO2)
    np.save(f"x_N2_{operation_step}",x_N2)
    np.save(f"Tw_ad_{operation_step}",Tw_ad)
    #======================================== Results - End
    
    result_data = result_calc(t,variable_vec)

    variable_vec = np.concatenate([P_ad[:,-1],T_ad[:,-1],y_CO2[:,-1],x_CO2[:,-1],x_N2[:,-1],Tw_ad[:,-1]])
    return variable_vec,result_data




def blowdown(variable_vec):
    global operation_step
    operation_step = 'blowdown'

    time_interval = np.array([0,t_blow])*v0/L

    #Solve_IVP ==============================================
    start = time.time()
    sol = scipy.integrate.solve_ivp(ode_system, time_interval, variable_vec, method=solve_ivp_method, rtol=1e-8, atol=1e-10)
    end = time.time()
    print(f"\nProcessing Time: {end - start:.2f} s")

    #Results ==============================================
    t = sol.t
    variable_vec = sol.y
    P_ad,T_ad,y_CO2,x_CO2,x_N2,Tw_ad = np.split(variable_vec, 6)

    #Velocity result in outlet ==============================================
    N_time_steps = len(t)
    v_ad = np.zeros(N_time_steps)
    for i in range(N_time_steps):
        P_ad_wall = np.zeros(N+1)
        P_ad_wall = WENO(P_ad[:,i],P_ad_wall)
        P_ad_wall[0] = P_ad[0,i] - v_ad_feed/dPdz_aux/2
        P_ad_wall[-1] = 1
        v_ad[i] = dPdz_aux*(P_ad_wall[-1]-P_ad_wall[-2])

    #Save ==============================================
    np.save(f"t_{operation_step}",t)
    np.save(f"P_ad_{operation_step}",P_ad)
    np.save(f"T_ad_{operation_step}",T_ad)
    np.save(f"v_ad_{operation_step}",v_ad)
    np.save(f"y_CO2_{operation_step}",y_CO2)
    np.save(f"x_CO2_{operation_step}",x_CO2)
    np.save(f"x_N2_{operation_step}",x_N2)
    np.save(f"Tw_ad_{operation_step}",Tw_ad)
    #======================================== Results - End
    
    result_data = result_calc(t,variable_vec)

    variable_vec = np.concatenate([P_ad[:,-1],T_ad[:,-1],y_CO2[:,-1],x_CO2[:,-1],x_N2[:,-1],Tw_ad[:,-1]])
    return variable_vec,result_data





def evacuation(variable_vec):
    global operation_step
    operation_step = 'evacuation'

    time_interval = np.array([0,t_evac])*v0/L

    #Solve_IVP ==============================================
    start = time.time()
    sol = scipy.integrate.solve_ivp(ode_system, time_interval, variable_vec, method=solve_ivp_method, rtol=1e-8, atol=1e-10)
    end = time.time()
    print(f"\nProcessing Time: {end - start:.2f} s")

    #Results ==============================================
    t = sol.t
    variable_vec = sol.y
    P_ad,T_ad,y_CO2,x_CO2,x_N2,Tw_ad = np.split(variable_vec, 6)

    #Velocity result in outlet ==============================================
    N_time_steps = len(t)
    v_ad = np.zeros(N_time_steps)
    for i in range(N_time_steps):
        P_ad_wall = np.zeros(N+1)
        P_ad_wall = WENO(P_ad[:,i],P_ad_wall)
        P_ad_wall[0] = P_ad[0,i] - v_ad_feed/dPdz_aux/2
        P_ad_wall[-1] = 1
        v_ad[i] = dPdz_aux*(P_ad_wall[-1]-P_ad_wall[-2])

    #Save ==============================================
    np.save(f"t_{operation_step}",t)
    np.save(f"P_ad_{operation_step}",P_ad)
    np.save(f"T_ad_{operation_step}",T_ad)
    np.save(f"v_ad_{operation_step}",v_ad)
    np.save(f"y_CO2_{operation_step}",y_CO2)
    np.save(f"x_CO2_{operation_step}",x_CO2)
    np.save(f"x_N2_{operation_step}",x_N2)
    np.save(f"Tw_ad_{operation_step}",Tw_ad)
    #======================================== Results - End
    
    result_data = result_calc(t,variable_vec)

    variable_vec = np.concatenate([P_ad[:,-1],T_ad[:,-1],y_CO2[:,-1],x_CO2[:,-1],x_N2[:,-1],Tw_ad[:,-1]])
    return variable_vec,result_data

def pressurization(variable_vec):
    global operation_step
    operation_step = 'pressurization'

    time_interval = np.array([0,t_press])*v0/L

    #Solve_IVP ==============================================
    start = time.time()
    sol = scipy.integrate.solve_ivp(ode_system, time_interval, variable_vec, method=solve_ivp_method, rtol=1e-8, atol=1e-10)
    end = time.time()
    print(f"\nProcessing Time: {end - start:.2f} s")

    #Results ==============================================
    t = sol.t
    variable_vec = sol.y
    P_ad,T_ad,y_CO2,x_CO2,x_N2,Tw_ad = np.split(variable_vec, 6)

    #Velocity result in outlet ==============================================
    N_time_steps = len(t)
    v_ad = np.zeros(N_time_steps)
    for i in range(N_time_steps):
        P_ad_wall = np.zeros(N+1)
        P_ad_wall = WENO(P_ad[:,i],P_ad_wall)
        P_ad_wall[0] = P_ad[0,i] - v_ad_feed/dPdz_aux/2
        P_ad_wall[-1] = 1
        v_ad[i] = dPdz_aux*(P_ad_wall[-1]-P_ad_wall[-2])

    #Save ==============================================
    np.save(f"t_{operation_step}",t)
    np.save(f"P_ad_{operation_step}",P_ad)
    np.save(f"T_ad_{operation_step}",T_ad)
    np.save(f"v_ad_{operation_step}",v_ad)
    np.save(f"y_CO2_{operation_step}",y_CO2)
    np.save(f"x_CO2_{operation_step}",x_CO2)
    np.save(f"x_N2_{operation_step}",x_N2)
    np.save(f"Tw_ad_{operation_step}",Tw_ad)
    #======================================== Results - End
    
    result_data = result_calc(t,variable_vec)

    variable_vec = np.concatenate([P_ad[:,-1],T_ad[:,-1],y_CO2[:,-1],x_CO2[:,-1],x_N2[:,-1],Tw_ad[:,-1]])
    return variable_vec,result_data

#Arrays ==============================================
t_global = np.zeros(1)
P_ad_global = np.zeros(N)
T_ad_global = np.zeros(N)
y_CO2_global = np.zeros(N)
x_CO2_global = np.zeros(N)
x_N2_global = np.zeros(N)
Tw_ad_global = np.zeros(N)

#Initial Condition ==============================================
# P_ad_global[:]  = P_L/P0
# T_ad_global[:]  = 1
# y_CO2_global[:] = 1e-6
# concentration   = np.array([1e-6,1-1e-6])*P_L/(R*T0)
# alpha, x_star   = solid_phase_calc(concentration,T0)
# x_CO2_global[:] = 0
# x_N2_global[:]  = x_star[1]
# Tw_ad_global[:] = 1
# variable_vec_initial = np.concatenate([P_ad_global,T_ad_global,y_CO2_global,x_CO2_global,x_N2_global,Tw_ad_global])
P_ad_global[:]  = np.load('P_ad_evacuation.npy')
T_ad_global[:]  = np.load('T_ad_evacuation.npy')
y_CO2_global[:] = np.load('y_CO2_evacuation.npy')
x_CO2_global[:] = np.load('x_CO2_evacuation.npy')
x_N2_global[:]  = np.load('x_N2_evacuation.npy')
Tw_ad_global[:] = np.load('Tw_ad_evacuation.npy')
variable_vec_initial = np.concatenate([P_ad_global,T_ad_global,y_CO2_global,x_CO2_global,x_N2_global,Tw_ad_global])
#======================================== Initial Condition - End

t_ads = 15
t_blow = 30
t_evac = 40
t_press = 15

steps_name = ['pressurization','adsorption','blowdown','evacuation']
steps_time = [t_press,t_ads,t_blow,t_evac]
steps_func = [pressurization,adsorption,blowdown,evacuation]

stop_criterium = 1000
tol = 1

n_steps = len(steps_time)

steps_feed = ['pressurization','adsorption']
steps_product = ['evacuation']

mole_CO2_in     = np.zeros(n_steps)
mole_CO2_out    = np.zeros(n_steps)
mole_CO2_acc    = np.zeros(n_steps)
mole_out        = np.zeros(n_steps)

heat_in         = np.zeros(n_steps)
heat_gen        = np.zeros(n_steps)
heat_lost       = np.zeros(n_steps)


def cycles(variable_vec):
    global stop_criterium,tol,t_step,t_actual
    global P_L, P_I


    P_L = 0.9e5
    P_I = 0.95e5
    n_cycles = 789
    stop_criterium = True
    CSS_indicator_list = np.array([])
    while stop_criterium:
        n_cycles += 1
        print(f'Cycle: {n_cycles}')
        t_actual = 0
        for i in range(n_steps):
            print(f'Operation:\t{steps_name[i]}')
            t_step = steps_time[i]

            variable_vec,result_data = steps_func[i](variable_vec)
            mole_CO2_in[i]  = result_data[0]
            mole_CO2_out[i] = result_data[1]
            mole_CO2_acc[i] = result_data[2]
            mole_out[i]     = result_data[4]
            
            heat_in[i]      = result_data[5]
            heat_gen[i]     = result_data[6]
            heat_lost[i]    = result_data[7]
            t_actual += steps_time[i]

        print(mole_CO2_in)
        print(mole_CO2_out)
        print(mole_out)
        
        mole_CO2_feed = 0
        for i in range(n_steps):
            if steps_name[i] in steps_feed:
                mole_CO2_feed += mole_CO2_in[i]
        
        mole_CO2_prod = 0
        mole_prod = 0
        for i in range(n_steps):
            if steps_name[i] in steps_product:
                mole_CO2_prod += mole_CO2_out[i]
                mole_prod += mole_out[i]
        
        purity = mole_CO2_prod/mole_prod
        recovery = mole_CO2_prod/mole_CO2_feed

        error_MB_cycle = np.abs(np.sum(mole_CO2_in)-np.sum(mole_CO2_out)-np.sum(mole_CO2_acc))/np.abs(np.sum(mole_CO2_in))
        error_EB_cycle = np.abs(np.sum(heat_in)+np.sum(heat_gen)-np.sum(heat_lost))/np.abs(np.sum(heat_in))

        print(f'\t Cycle:  {n_cycles:3d}\tPurity:  {purity*100:4.1f}%\tRecovery:  {recovery*100:4.1f}%\tError(MB):  {error_MB_cycle*100:4.1f}%\tError(EB):  {error_EB_cycle*100:4.1f}%')#,end='\r')


        CSS_indicator = 1 - np.sum(mole_CO2_out)/np.sum(mole_CO2_in)
        print(f'{CSS_indicator*100:.2f}')
        if n_cycles > 50:
            CSS_diff = (CSS_indicator_list[-6:-1] - CSS_indicator)/CSS_indicator
            stop_criterium = False
            for i in range(5):
                if CSS_diff[i] > 0.005:
                    stop_criterium = True
                    break
            if stop_criterium == False:
                print('\tConvergence of CSS indicator!')

        CSS_indicator_list = np.append(CSS_indicator_list,CSS_indicator)
        np.save('CSS_indicator_list.npy',CSS_indicator_list)

        if n_cycles > 2000:
            stop_criterium = False
            print('\tMaximum number of cycles reached')

        n_stabilization_cycles = 400
        if n_cycles > n_stabilization_cycles:
            if P_L > P_L_final:
                P_L -= 0.1/n_stabilization_cycles*P0
            else:
                P_L = P_L_final
            if P_I > P_I_final:
                P_I -= 0.1/n_stabilization_cycles*P0
            else:
                P_I = P_I_final

    return

cycles(variable_vec_initial)