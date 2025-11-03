#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:51:28 2023

@author: mikinakajima
s"""

import matplotlib.pyplot as plt
#from matplotlib import font_manager
import numpy as np
from scipy import interpolate
import sys
from scipy.optimize import curve_fit

plt.rcParams['font.family'] = 'Helvetica'



filenames = ['MgFeO_shock_data/39878_Us_T.txt', 'MgFeO_shock_data/39874_Us_T.txt', 'MgFeO_shock_data/38691_Us_T.txt', 'MgFeO_shock_data/38692_Us_T.txt', 'MgFeO_shock_data/39879_Us_T.txt', 'MgFeO_shock_data/38694_Us_T.txt',    'MgFeO_shock_data/38693_Us_T.txt',  'MgFeO_shock_data/39877_Us_T.txt']
colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)','Mg$_{0.95}$Fe$_{0.05}$O (38694)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)',  'Mg$_{0.95}$Fe$_{0.05}$O (39877)'  ]


                                                            
# Reading DFT-MD data

file_DFT = 'Hugoniot/mgo_dft_data.txt'


                                                              
 
def open_files(filename, ax,labelname,color):

    data = np.loadtxt(filename)

    xx = data[:, 0] # shock velocity
    yy = data[:, 3] # Temperature

    x_fit = np.linspace(min(xx),max(xx),100)
    model = np.poly1d(np.polyfit(xx, yy, 2))
    error = np.sqrt(np.sum((model(xx) - yy)**2/len(xx)))


    ax.scatter(xx,yy, s=4, label=labelname,c=color,marker = 'o', alpha=0.3, edgecolors='none')
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



fig, ax1 = plt.subplots()


x2 = np.arange(16, 29)

ax1.set_xticks(x2)

#for i in range(0,len(filenames)):
for i in range(0,len(filenames)):
    open_files(filenames[i], ax1,labelnames[i],colors[i])
    

#combining MgFeO dataset
data = np.vstack([np.loadtxt(f) for f in filenames])
xx_combined = data[:, 0] # shock velocity
yy_combined = data[:, 3] # Temperature


x_fit = np.linspace(min(xx_combined),max(xx_combined),100)
model0, cov_matrix = np.polyfit(xx_combined, yy_combined, 2, cov=True)
model = np.poly1d(model0)
error = np.sqrt(np.sum((model(xx_combined) - yy_combined)**2/len(xx_combined)))

coefficients = model.coefficients
errors = np.sqrt(np.diag(cov_matrix))
print(coefficients)
print(errors)
np.savetxt("T_fit_coefficients.txt", coefficients)



#ax1.scatter(xx,yy, s=4, label=labelname,facecolor = color)
ax1.plot(x_fit, model(x_fit),color='grey',alpha=0.2,zorder=1,label='Fit',linewidth='0.5')
ax1.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor='grey', alpha=0.1,zorder=1)

ax1.errorbar(V_Mc, T_Mc, xerr = V_Mc_error, yerr=T_Mc_error, fmt='none', color='tan', label='MgO, M (2019)', elinewidth=0.9)

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


ax1.set_xlabel('Shock Velocity (km/s)', fontsize=17)
ax1.set_ylabel('Temperature ($10^3$ K)', fontsize=17)
ax2.set_xlabel('Pressure (GPa)', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)


#ax1.legend(frameon=False,fontsize=8)
fig.savefig('Fig_Us_T.pdf',dpi=1000,  bbox_inches="tight")
plt.show()
plt.close()
