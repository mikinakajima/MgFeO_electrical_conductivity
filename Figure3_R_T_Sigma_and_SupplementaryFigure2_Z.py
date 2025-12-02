#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating Figure 3 of Nakajima et al. , entitled "Electrical conductivities of (Mg,Fe)O at extreme pressures and implications for planetary magma oceans"

"""

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42   # Embed text as TrueType (editable)
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Helvetica'


import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.constants import pi, G, hbar, m_e, m_u, k, e, N_A



fig, axes = plt.subplots(3, 1, figsize=(7.2, 10.2), sharex=True)
fig.subplots_adjust(hspace=0.3)

# ---------- Figure 3 PANEL (a): Reflectivity ----------

axes[0].tick_params(axis='x', which='both', labelbottom=True)  
axes[1].tick_params(axis='x', which='both', labelbottom=True)


filenames = ['MgFeO_shock_data/39878_Rfit.txt', 'MgFeO_shock_data/39874_Rfit.txt', 'MgFeO_shock_data/38691_Rfit.txt', 'MgFeO_shock_data/38692_Rfit.txt', 'MgFeO_shock_data/39879_Rfit.txt', 'MgFeO_shock_data/38694_Rfit.txt', 'MgFeO_shock_data/39882_Rfit.txt',   'MgFeO_shock_data/38693_Rfit.txt',  'MgFeO_shock_data/39877_Rfit.txt']

colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue', 'navy' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','(Mg$_{0.98}$,Fe$_{0.02}$)O (38692)', '(Mg$_{0.98}$,Fe$_{0.02}$)O (39879)','(Mg$_{0.95}$,Fe$_{0.05}$)O (38694)', '(Mg$_{0.95}$,Fe$_{0.05}$)O (39882)', '(Mg$_{0.95}$,Fe$_{0.05}$)O (38693)',  '(Mg$_{0.95}$,Fe$_{0.05}$)O (39877)']
                                       
 
   # Soubiran & Militzer (2018)
SMdata_P = [556.425, 631.8, 691.3, 754.38]  #in GPa
SMdata_R = [1.656, 3.017, 4.568, 6.63]      #reflectivity    
#SMdata_T_P = [554.80,628.583, 692.179, 754.495]
SMdata_T = np.array([11.959,14.046, 15.983, 18.032])*1000.0




                                                       
 
def open_files(filename, ax,labelname,color):
    
    data = np.loadtxt(filename)
    xx = data[:, 1] # shock velocity from our experiments
    yy = data[:, 0] # reflectivity from our experiments

    x_fit = np.linspace(xx.min(), xx.max(), 100)
    
    
    #modifications to the models to the two dataset so that the fit increases as the pressure increases
    if filename == 'MgFeO_shock_data/39879_Rfit.txt' or filename == 'MgFeO_shock_data/39874_Rfit.txt':
        model = np.poly1d(np.polyfit(xx, yy, 3))
    else:
        model = np.poly1d(np.polyfit(xx, yy, 2))
        
    error = np.sqrt(np.sum((model(xx) - yy)**2/len(xx)))
    ax.scatter(xx,yy, s=4,label=labelname, c=color, marker = 'o', alpha=0.3, edgecolors='none', zorder=3)
    ax.plot(x_fit, model(x_fit),color=color,zorder=3)
    ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2,zorder=3)



x = np.arange(16, 29)


rho0, a, b = np.loadtxt("Up-Us_coefficients.txt")
b = b *1000.0


Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments

SMdata_Us = []
f = interpolate.interp1d(Press, x)
for i in range(0,len(SMdata_P)):
    SMdata_Us.append(f(SMdata_P[i]))

#fig, ax1 = plt.subplots()
ax1 = axes[0]
ax1.set_xticks(x)

for i in range(0,9):
    open_files(filenames[i], ax1,labelnames[i],colors[i])

# Reading DFT-MD data
file_DFT = 'Hugoniot/mgo_dft_data.txt'
data = np.loadtxt(file_DFT, skiprows=1)
#shock velocity, reflectivity, reflectivity error from DFT-MD
V_DFT, Ref_DFT, Ref_DFT_error = data[:, 4], data[:, 7]*100.0, data[:, 8]*100.0 


data = np.loadtxt('Hugoniot/McCoy_R.txt', skiprows=1)
#shock velocity, shock velocity error, reflectivity, reflectivity error from McCoy et al 2019
V_Mc, V_Mc_error, Ref_Mc, Ref_Mc_error = data[:, 1], data[:, 2], data[:, 3] * 100.0, data[:, 4] * 100.0


def polynomials(x,a,b,c,d):
    return a * x**3.0 + b*x**2.0 + c*x +d



ax1.errorbar(V_DFT[0], Ref_DFT[0], yerr=Ref_DFT_error[0], fmt='x', color='blue', label='MgO, DFT-MD (B2)', zorder=3)
ax1.errorbar(V_DFT[1:5], Ref_DFT[1:5], yerr=Ref_DFT_error[1:5], fmt='x', color='red', label='MgO, DFT-MD (Liquid)', zorder=3)
ax1.scatter(SMdata_Us, SMdata_R, label='MgO, DFT-MD (Soubiran & Militzer 2018)', zorder=1)
ax1.errorbar(V_Mc, Ref_Mc, xerr = V_Mc_error, yerr=Ref_Mc_error, fmt='none', color='tan', label='MgO (McCoy et al. 2019)', elinewidth=0.9, zorder=1)


ax1.set_xlim([15.5, 28])
ylimit = [-1,25]
ax1.set_ylim([ylimit[0],ylimit[1] ])
ax1.text(26.5, ylimit[0]+0.1*(ylimit[1]-ylimit[0]),'(a)' ,fontsize=20)


ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1.get_xticks())
Press_int = [int(value) for value in Press]
ax2.set_xticklabels(Press_int)


# Set custom ticks on the top x-axis (corresponding pressure values)
custom_press_ticks = [400, 600, 800, 1000, 1200, 1400, 1600]

# Interpolate corresponding x-axis values for the custom pressure ticks
from scipy.interpolate import interp1d
pressure_to_velocity = interp1d(Press, x, bounds_error=False, fill_value="extrapolate")
velocity_for_press_ticks = pressure_to_velocity(custom_press_ticks)

# Set ticks and labels on the top axis
ax2.set_xticks(velocity_for_press_ticks)
ax2.set_xticklabels(custom_press_ticks)


for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)


#ax1.set_xlabel('Shock Velocity (km/s)', fontsize=17)
ax1.set_ylabel('Reflectivity, $R$ (%)', fontsize=17)
ax2.set_xlabel('Pressure (GPa)', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)

handles, labels = ax1.get_legend_handles_labels()

# --- split into the two groups ---
exp_handles = handles[:9]    # our data 
exp_labels  = labels[:9]

lit_handles = handles[9:13]  # 9–12 literature data
lit_labels  = labels[9:13]

# desired literature order: B2 → Liquid → Soubiran → McCoy
lit_order = [1, 2, 0, 3]
lit_handles = [lit_handles[i] for i in lit_order]
lit_labels  = [lit_labels[i]  for i in lit_order]

# --- combine and make legend ---
all_handles = exp_handles + lit_handles
all_labels  = exp_labels  + lit_labels

ax1.legend(all_handles, all_labels, frameon=False, fontsize=6)



# ---------- Figure 3 PANEL (b): Temperature ----------



filenames = ['MgFeO_shock_data/39878_Us_T.txt', 'MgFeO_shock_data/39874_Us_T.txt', 'MgFeO_shock_data/38691_Us_T.txt', 'MgFeO_shock_data/38692_Us_T.txt', 'MgFeO_shock_data/39879_Us_T.txt', 'MgFeO_shock_data/38694_Us_T.txt',    'MgFeO_shock_data/38693_Us_T.txt',  'MgFeO_shock_data/39877_Us_T.txt']
colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)','Mg$_{0.95}$Fe$_{0.05}$O (38694)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)',  'Mg$_{0.95}$Fe$_{0.05}$O (39877)'  ]


                                                            
# Reading DFT-MD data

file_DFT = 'Hugoniot/mgo_dft_data.txt'


                                                              
 
def open_files(filename, ax,labelname,color):

    data = np.loadtxt(filename)

    xx = data[:, 0] # shock velocity
    yy = data[:, 1] # Temperature

    x_fit = np.linspace(min(xx),max(xx),100)
    model = np.poly1d(np.polyfit(xx, yy, 2))
    error = np.sqrt(np.sum((model(xx) - yy)**2/len(xx)))


    ax.scatter(xx,yy, s=4, label=labelname,c=color,marker = 'o', alpha=0.3, edgecolors='none',zorder=2)
    ax.plot(x_fit, model(x_fit),color=color,zorder=2)
    ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2,zorder=2)

x = np.arange(15, 31)



rho0, a, b = np.loadtxt("Up-Us_coefficients.txt")
b = b *1000.0


Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments

f = interpolate.interp1d(Press, x)
      
# liquid-B2 MgO phase curve
data = np.loadtxt('phase/phase.txt', delimiter=',')
mask = data[:, 0] > 344.0 # mask low pressures ranges so that interpolation is easier

Phase_P = data[mask, 0]
Phase_T = data[mask, 1] * 1000.0
Phase_Us = f(Phase_P)
            
# B1-B2 phase curve            
phase_file = 'phase/phase_B1.txt'
data = np.loadtxt(phase_file, delimiter=',')
mask = data[:, 0] > 344.0

# Extract filtered columns
Phase_P_B1 = data[mask, 0]
Phase_T_B1 = data[mask, 1] * 1000.0
Phase_Us_B1 = f(Phase_P_B1)
            

# Function to interpolate temperature data
def interpolate_temperature(pressure_values, temperature_values):
    # Convert lists to numpy arrays for easier manipulation
    pressure_values = np.array(pressure_values)
    temperature_values = np.array(temperature_values)
    
    # Fit a polynomial of degree 3 (adjust degree as needed) to the data
    poly_degree = 3
    coeffs = np.polyfit(pressure_values, temperature_values, poly_degree)
    
    # Create a polynomial function from the coefficients
    poly_func = np.poly1d(coeffs)
    
    # Generate a smooth curve using the polynomial function
    pressure_smooth = np.linspace(min(pressure_values), max(pressure_values), 100)
    temperature_smooth = poly_func(pressure_smooth)
    
    return pressure_smooth, temperature_smooth



# Interpolate temperature data
B1_Us_smooth, B1_temperature_smooth = interpolate_temperature(Phase_Us_B1, Phase_T_B1)
Us_smooth, T_smooth = interpolate_temperature(Phase_Us, Phase_T)


#McCoy et al 2019
data = np.loadtxt('Hugoniot/McCoy_R.txt', skiprows=1)

# velocity, velocity error, tempearture, temperature error from McCoy et al 2019
V_Mc        = data[:, 1]
V_Mc_error  = data[:, 2]
T_Mc        = data[:, 5]*1000.0
T_Mc_error  = data[:, 6]*1000.0




ax1 = axes[1]


x2 = np.arange(16, 29)

ax1.set_xticks(x2)

for i in range(0,len(filenames)):
    open_files(filenames[i], ax1,labelnames[i],colors[i])
    

#combining MgFeO dataset
data = np.vstack([np.loadtxt(f) for f in filenames])
xx_combined = data[:, 0] # shock velocity
yy_combined = data[:, 1] # Temperature


x_fit = np.linspace(min(xx_combined),max(xx_combined),100)
model0, cov_matrix = np.polyfit(xx_combined, yy_combined, 2, cov=True)
model = np.poly1d(model0)
error = np.sqrt(np.sum((model(xx_combined) - yy_combined)**2/len(xx_combined)))

coefficients = model.coefficients
errors = np.sqrt(np.diag(cov_matrix))
print(coefficients)
print(errors)
np.savetxt("T_fit_coefficients.txt", coefficients)



ax1.plot(x_fit, model(x_fit),color='grey',alpha=0.2,zorder=1,label='Fit',linewidth='0.5')
ax1.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor='grey', alpha=0.1,zorder=1)

ax1.errorbar(V_Mc, T_Mc, xerr = V_Mc_error, yerr=T_Mc_error, fmt='none', color='tan', label='MgO, M (2019)', elinewidth=0.9, zorder = 1)

ax1.scatter(SMdata_Us, SMdata_T, label='MgO, DFT-MD (Soubiran & Militzer 2018)', zorder=1)

#producing output file that has P-T information
file_P_T = 'P_T.txt'
f = open(file_P_T, 'w')


Press_combined = rho0 * xx_combined * 1000 * (xx_combined * 1000.0-b)/a * 1e-9 #Us-P relationship based on our experiments
for i in range(0,len(xx_combined)):
    f.write(f"{Press_combined[i]} {yy_combined[i]}\n")
f.close()


# Load DFT-MD data, skipping the header line
data = np.loadtxt(file_DFT, skiprows=1)

# Extract columns
V_DFT = data[:, 4]
Temp_DFT = data[:, 2]
        

ax1.errorbar(V_DFT[0], Temp_DFT[0], fmt='x', color='blue', label='MgO, DFT(B2)')
ax1.errorbar(V_DFT[1:5], Temp_DFT[1:5], fmt='x', color='red', label='MgO, DFT(Liquid)')

ax1.plot(Us_smooth, T_smooth,linewidth=0.5,color='steelblue')
ax1.fill_between(Phase_Us,Phase_T, 60000, color='seashell',zorder=0)
ax1.fill_between(Phase_Us,0, Phase_T,color='aliceblue',zorder=0)
ax1.plot(B1_Us_smooth, B1_temperature_smooth,linewidth=0.5,color='gold')
ax1.fill_between(B1_Us_smooth,0, B1_temperature_smooth,color='cornsilk',zorder=0)
ax1.text(15.6, 2000, 'B1 (MgO)',color='gold',fontsize=10)
ax1.text(19, 5000, 'B2 (MgO)',color='steelblue',fontsize=10)
ax1.text(24, 20000, 'Liquid (MgO)',color='lightcoral',fontsize=10)

for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)

ax1.set_xlim([15.5, 28])
ylimit = [0,60000]
ax1.set_ylim([ylimit[0],ylimit[1] ])
ax1.text(26.5, ylimit[0]+0.1*(ylimit[1]-ylimit[0]),'(b)' ,fontsize=20)

y = np.array([0, 1e4, 2e4, 3e4, 4e4, 5e4, 6e4 ])
y_tick_labels= [f'$0$', f'$10$',f'$20$', f'$30$' , f'$40$', f'$50$', f'$60$']

ax1.set_yticks(y)
ax1.set_yticklabels(y_tick_labels)


ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1.get_xticks())


Press = rho0 * x2 * 1000 * (x2 * 1000-b)/a * 1e-9
Press_int = [int(value) for value in Press]
ax2.set_xticklabels(Press_int)

# Set custom ticks on the top x-axis (corresponding pressure values)
custom_press_ticks = [400, 600, 800, 1000, 1200, 1400, 1600]

# Interpolate corresponding x-axis values for the custom pressure ticks
from scipy.interpolate import interp1d
pressure_to_velocity = interp1d(Press, x2, bounds_error=False, fill_value="extrapolate")
velocity_for_press_ticks = pressure_to_velocity(custom_press_ticks)

# Set ticks and labels on the top axis
ax2.set_xticks(velocity_for_press_ticks)
ax2.set_xticklabels(custom_press_ticks)

ax1.set_ylabel('Temperature ($10^3$ K)', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)





# ---------- Figure 3 PANEL (C): DC Conductivity ----------



plt.rcParams['font.family'] = 'Helvetica'
filenames = ['MgFeO_shock_data/39878_Rfit.txt', 'MgFeO_shock_data/39874_Rfit.txt', 'MgFeO_shock_data/38691_Rfit.txt', 'MgFeO_shock_data/38692_Rfit.txt', 'MgFeO_shock_data/39879_Rfit.txt', 'MgFeO_shock_data/38694_Rfit.txt', 'MgFeO_shock_data/39882_Rfit.txt',   'MgFeO_shock_data/38693_Rfit.txt',  'MgFeO_shock_data/39877_Rfit.txt']

colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue', 'navy' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)','Mg$_{0.95}$Fe$_{0.05}$O (38694)', 'Mg$_{0.95}$Fe$_{0.05}$O (39882)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)',  'Mg$_{0.95}$Fe$_{0.05}$O (39877)'  ]



file_DFT = 'Hugoniot/mgo_dft_data.txt' #DFT data
# generating a file that includes Us-DC conductivity file
file_Us_sigma = 'Us_sigma.txt'

n0 = 1.743 #initial refractive index of MgO from McWilliams 2012
e_0 = 8.854e-12 #Vacuum permittivity
e_charge=1.602e-19 #electron charge, C, A/s
molmass = 40.0 * 1e-3  #kg/mol, MgO, assuming no dissociation

rho0, a, b = np.loadtxt("Up-Us_coefficients.txt")
b = b *1000.0
a_temp, b_temp, c_temp = np.loadtxt("T_fit_coefficients.txt")

f = open(file_Us_sigma, 'w')
# writing a file for the shock velocity vs electrical conductivity

def open_files(filename, labelname, color, ax, axZ):

        
    # ----- Load data -----
# assuming 2-column file: RR(%) and Us(km/s)
    data = np.loadtxt(filename)
    RR = data[:, 0] * 0.01       # reflectivity, dimensionless
    Us = data[:, 1]               # shock velocity, km/s

# ----- Derived quantities -----
    up   = (Us * 1000 - b) / a                             # particle velocity [m/s]
    rho  = rho0 * Us * 1000 / (Us * 1000 - up)             # shocked density [kg/m^3]
    temp = a_temp * Us**2 + b_temp * Us + c_temp            # temperature [K or as defined]
    pressure = rho0 * Us * 1000 * up * 1e-9                 # pressure [GPa]
        

    sigma = []
    Z_store = []
    
    
    for i in range(0,len(rho)):
        Z_ini = 0.0001
        Z = Z_ini # ionization, free parameter  
        r_nb = 0.0
        z_l = 0.0
        z_h = 0.0
        
        
        if RR[i]<=0.0:
            continue
        else:
    
            while (RR[i]-r_nb) > 0:
                r_nb, sigma0, omegatau0 = calculate_r_nb(rho[i],Z,temp[i])
                Z = Z + 0.01
            z_l = Z - 0.02
            if z_l <0:
                z_l = Z_ini
            z_h = Z - 0.01
        
            Z = (z_l + z_h)*0.5
            r_nb, sigma0, omegatau0 = calculate_r_nb(rho[i],Z,temp[i])
        
        
            while abs(RR[i]-r_nb)/RR[i]>0.01:
                r_nb, sigma0, omegatau0 = calculate_r_nb(rho[i],Z,temp[i])    
        
                if RR[i]-r_nb>0:
                    z_l = Z
                    z_h = z_h
                else:
                    z_l = z_l
                    z_h = Z
        
                Z = (z_l + z_h)*0.5
                r_nb, sigma0, omegatau0 = calculate_r_nb(rho[i],Z,temp[i])


                r_nb_l, sigma1, omegatau1 = calculate_r_nb(rho[i],z_l,temp[i])
                r_nb_h, sigma1,ometagau1 = calculate_r_nb(rho[i],z_h,temp[i])
            
            
        sigma.append(sigma0)
        Z_store.append(Z)
    
        
    Us_filtered = []
    sigma_filtered = []
    pressure_filtered = []
    Z_filtered = []

    for i in range(len(Us)):
        if Us[i] > 18:
            Us_filtered.append(Us[i])
            sigma_filtered.append(sigma[i])
            pressure_filtered.append(pressure[i])
            Z_filtered.append(Z_store[i])            
            
    ax.scatter(Us_filtered,sigma_filtered, s=4, c=color, marker = 'o', alpha=0.3, edgecolors='none') 
    axZ.scatter(Us_filtered,Z_filtered, s=4, c=color, marker = 'o', alpha=0.3, edgecolors='none') 

    
    for i in range(0,len(Us_filtered)):
        f.write(f"{Us_filtered[i]} {sigma_filtered[i]} {pressure_filtered[i]}\n")
        

    if len(Us_filtered)>18.0: # Figure 3
        x_fit = np.linspace(min(Us_filtered),max(Us_filtered),100)
        model = np.poly1d(np.polyfit(Us_filtered, sigma_filtered, 2))
        error = np.sqrt(np.sum((model(Us_filtered) - sigma_filtered)**2/len(Us_filtered)))
        ax.plot(x_fit, model(x_fit),color=color)
        ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2)



    if len(Us_filtered)>18.0: # Supplementary Figure 2 (Z values)
        x_fit = np.linspace(min(Us_filtered),max(Us_filtered),100)
        model = np.poly1d(np.polyfit(Us_filtered, Z_filtered, 2))
        error = np.sqrt(np.sum((model(Us_filtered) - Z_filtered)**2/len(Us_filtered)))
        axZ.plot(x_fit, model(x_fit),color=color)
        axZ.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2)

 

def calculate_r_nb(rho,Z, temp):


    nb =  1.598 # from Figure S3 McWilliams et al 2012

    ni = (rho*N_A)/molmass #ion density
    ne = Z*ni #electron density 

    v_fermi = (hbar/m_e) * (3 * pi**2 *ne)**(1.0/3.0) #fermi velocity 
    v_thermal = np.sqrt((2.0*k*temp)/m_e) #thermal velocity
    vel = [v_fermi, v_thermal] # use the larger velocity between fermi and thermal velocities
    l =  2.0 * (3.0/(4.0*np.pi*ni))**(1.0/3.0)

    omega = 2.0 * np.pi * 2.99792e8/(532e-9) 
    omega_p = np.sqrt((ne * e_charge**2)/(e_0*m_e)) #plasma frequency
    scatter_tau = l/max(vel)
    
    omegatau=scatter_tau * omega
    
    
    sigma0 = (ne * e_charge **2 * scatter_tau)/m_e #DC conductivity eq. 43 from Millot 2015


    n_mod =np.sqrt(nb**2 - (omega_p**2) / (omega**2 * (1 + 1j/(omega*scatter_tau))))

    
    r_nb = (np.abs(n0-n_mod)/np.abs(n0+n_mod))**2.0


    return r_nb, sigma0, omegatau

ax1 = axes[2]
figZ, axZ = plt.subplots()

for i in range(0,len(filenames)):
    open_files(filenames[i], labelnames[i],colors[i],ax1,axZ)
 
    

x = np.arange(16,29)
ax1.set_xticks(x)

y = np.array([0, 5e4, 10e4, 15e4, 20e4, 25e4, 30e4 ])
y_tick_labels= [f'$0$', f'$5$',f'$10$', f'$15$' , f'$20$', f'$25$', f'$30$']

ax1.set_yticklabels(y_tick_labels)
ax1.set_yticks(y)


f.close()

ax1.set_xlim([15.5, 28])


ylimit = [-10000,250000]
ax1.set_ylim([ylimit[0],ylimit[1]])
ax1.text(26.5, ylimit[0]+0.1*(ylimit[1]-ylimit[0]),'(c)' ,fontsize=20)

Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments
Temp =  a_temp * x**2.0 + b_temp*x + c_temp



# Load data from file, skipping the header
data = np.loadtxt(file_DFT, skiprows=1)

V_DFT = data[:, 4]
Sigma_DFT = data[:, 5]
Sigma_DFT_error = data[:, 6]


ax1.errorbar(V_DFT[0], Sigma_DFT[0], yerr=Sigma_DFT_error[0], fmt='x', color='blue')#, label='DFT, MgO(B2)')
ax1.errorbar(V_DFT[1:5], Sigma_DFT[1:5], yerr=Sigma_DFT_error[1:5], fmt='x', color='red')#, label='DFT, MgO(Liquid)')



ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1.get_xticks())
Press_int = [int(value) for value in Press]
ax2.set_xticklabels(Press_int)


# Set custom ticks on the top x-axis (corresponding pressure values)
custom_press_ticks = [400, 600, 800, 1000, 1200, 1400, 1600]

# Interpolate corresponding x-axis values for the custom pressure ticks

pressure_to_velocity = interp1d(Press, x, bounds_error=False, fill_value="extrapolate")
velocity_for_press_ticks = pressure_to_velocity(custom_press_ticks)

# Set ticks and labels on the top axis
ax2.set_xticks(velocity_for_press_ticks)
ax2.set_xticklabels(custom_press_ticks)





# strixrude + 2020 comparison
sigma00 = 1.99e5

Delta_E = 75.94e3
Delta_V = -0.061 * 1e-6
R = 8.314 # gas constant
P = Press * 1e9

sigma_stix=np.zeros(len(x))
for i in range(0,len(x)):
    sigma_stix[i] = sigma00 * np.exp(-(Delta_E+P[i]*Delta_V)/(R*Temp[i]))

ax1.plot(x, sigma_stix*10**(-4), label='MgO ', color = 'tan', linestyle='-')
ax1.plot(x, sigma_stix*10**(-3.2), label='Mg$_{0.95}$Fe$_{0.05}$O ', color = 'navy', linestyle='dashed')
ax1.plot(x, sigma_stix, label='Mg$_{0.75}$Fe$_{0.25}$O ', color = 'grey', linestyle='dashed')
ax1.text(15.7, 23.2*1e4, 'DFT-MD (Holmström et al. 2018)', fontsize=6)

ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)




data = np.loadtxt(file_Us_sigma)
xx_mod = data[:, 0]  # Shock velocity
yy_mod = data[:, 1]  # Reflectivity



x_fit = np.linspace(min(xx_mod),max(xx_mod),100)
model0, cov_matrix = np.polyfit(xx_mod, yy_mod, 2, cov=True) #np.poly1d(np.polyfit(xx_mod, yy_mod, 2 ))
model = np.poly1d(model0)
error = np.sqrt(np.sum((model(xx_mod) - yy_mod)**2/len(xx_mod)))

coefficients = model.coefficients
errors = np.sqrt(np.diag(cov_matrix))


x_fit2 = np.linspace(min(xx_mod),max(xx_mod),100)

print("Coefficients:", coefficients)
print("errors:", errors)


for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)


ax1.set_xlabel('Shock Velocity (km/s)', fontsize=17)
ax1.set_ylabel('DC Conductivity, $\sigma$ ($10^4$S/m)', fontsize=17)
ax1.legend(frameon=False,fontsize=6, loc = "upper left", bbox_to_anchor=(0.01, 0.95))




fig.savefig('Figure3.pdf',dpi=1000,  bbox_inches="tight")




# ---------- Supplementary Figure 2: Z values ----------


    
x = np.arange(16,29)
axZ.set_xticks(x)

y = np.array([0, 0.1, 0.2,0.3, 0.4])
axZ.set_yticks(y)

axZ.set_xlim([15.5, 28])

ylimit = [0,0.35]
axZ.set_ylim([ylimit[0],ylimit[1]])
Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments
Temp =  a_temp * x**2.0 + b_temp*x + c_temp


axZ2 = axZ.twiny()
axZ2.set_xlim(axZ.get_xlim())
axZ2.set_xticks(axZ.get_xticks())
Press_int = [int(value) for value in Press]
axZ2.set_xticklabels(Press_int)


# Set custom ticks on the top x-axis (corresponding pressure values)
custom_press_ticks = [400, 600, 800, 1000, 1200, 1400, 1600]

# Interpolate corresponding x-axis values for the custom pressure ticks
from scipy.interpolate import interp1d
pressure_to_velocity = interp1d(Press, x, bounds_error=False, fill_value="extrapolate")
velocity_for_press_ticks = pressure_to_velocity(custom_press_ticks)

axZ2.set_xticks(velocity_for_press_ticks)
axZ2.set_xticklabels(custom_press_ticks)


for i in ['top', 'bottom','left','right']:
    axZ.spines[i].set_linewidth(1.5)


axZ.tick_params(axis='both', which='major', labelsize=11)
axZ2.tick_params(axis='both', which='major', labelsize=11)


axZ.set_xlabel('Shock Velocity (km/s)', fontsize=16)
axZ.set_ylabel('Z values', fontsize=16)
axZ2.set_xlabel('Pressure (GPa)', fontsize=16)

for i in ['top', 'bottom','left','right']:
    axZ.spines[i].set_linewidth(1.5)


axZ.legend(frameon=False,fontsize=8, loc = "upper left", bbox_to_anchor=(0.01, 0.95))
figZ.savefig('SupplementaryFigure2.pdf',dpi=1000)

plt.show()
plt.close()






