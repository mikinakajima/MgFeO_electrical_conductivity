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
from scipy.optimize import curve_fit

plt.rcParams['font.family'] = 'Helvetica'


filenames = ['MgFeO_shock_data/39878_Rfit.txt', 'MgFeO_shock_data/39874_Rfit.txt', 'MgFeO_shock_data/38691_Rfit.txt', 'MgFeO_shock_data/38692_Rfit.txt', 'MgFeO_shock_data/39879_Rfit.txt', 'MgFeO_shock_data/38694_Rfit.txt', 'MgFeO_shock_data/39882_Rfit.txt',   'MgFeO_shock_data/38693_Rfit.txt',  'MgFeO_shock_data/39877_Rfit.txt']

colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue', 'navy' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','(Mg$_{0.98}$,Fe$_{0.02}$)O (38692)', '(Mg$_{0.98}$,Fe$_{0.02}$)O (39879)','(Mg$_{0.95}$,Fe$_{0.05}$)O (38694)', '(Mg$_{0.95}$,Fe$_{0.05}$)O (39882)', '(Mg$_{0.95}$,Fe$_{0.05}$)O (38693)',  '(Mg$_{0.95}$,Fe$_{0.05}$)O (39877)']
                                       
 
   # Soubiran & Militzer (2018)
SMdata_P = [556.425, 631.8, 691.3, 754.38]  #in GPa
SMdata_R = [1.656, 3.017, 4.568, 6.63]      #reflectivity                                                #Burkhardt

                                                              
 
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
    ax.scatter(xx,yy, s=4,label=labelname, c=color, marker = 'o', alpha=0.3, edgecolors='none')
    ax.plot(x_fit, model(x_fit),color=color)
    ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2)



x = np.arange(16, 29)


rho0, a, b = np.loadtxt("Up-Us_coefficients.txt")
b = b *1000.0


Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments

SMdata_Us = []
f = interpolate.interp1d(Press, x)
for i in range(0,len(SMdata_P)):
    SMdata_Us.append(f(SMdata_P[i]))

fig, ax1 = plt.subplots()
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



ax1.errorbar(V_DFT[0], Ref_DFT[0], yerr=Ref_DFT_error[0], fmt='x', color='blue', label='MgO, DFT-MD (B2)')
ax1.errorbar(V_DFT[1:5], Ref_DFT[1:5], yerr=Ref_DFT_error[1:5], fmt='x', color='red', label='MgO, DFT-MD (Liquid)')
ax1.scatter(SMdata_Us, SMdata_R, label='MgO, DFT-MD (Soubiran & Militzer 2018)')

ax1.errorbar(V_Mc, Ref_Mc, xerr = V_Mc_error, yerr=Ref_Mc_error, fmt='none', color='tan', label='MgO (McCoy et al. 2019)', elinewidth=0.9)


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


ax1.set_xlabel('Shock Velocity (km/s)', fontsize=17)
ax1.set_ylabel('Reflectivity, $R$ (%)', fontsize=17)
ax2.set_xlabel('Pressure (GPa)', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)

ax1.legend(frameon=False,fontsize=8)
fig.savefig('Fig_Rfit.pdf',dpi=1000, bbox_inches="tight")
plt.show()
plt.close()
