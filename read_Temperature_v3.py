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

# to add more data, just add a filename, label name, and your favorite color.
#filenames = ['../omegaMgO_38691/38691_Us_T.txt','../omegaMgFeO_38692/38692_Us_T.txt', '../omegaMgFeO_38693/38693_Us_T.txt','../omegaMgFeO_38694/38694_Us_T.txt','../omegaMgO_39874/39874_Us_T.txt','../omegaMgFeO_39877/39877_Us_T.txt','../omegaMgO_39878/39878_Us_T.txt', '../omegaMgFeO_39879/39879_Us_T.txt']  

#labelnames = ['MgO', 'Mg$_{0.98}$Fe$_{0.02}$O', 'Mg$_{0.95}$Fe$_{0.05}$O', 'Mg$_{0.95}$Fe$_{0.05}$O', 'MgO', 'Mg$_{0.95}$Fe$_{0.05}$O','MgO','Mg$_{0.95}$Fe$_{0.05}$O']#, 'MgO', 'Mg$_{0.98}$Fe$_{0.02}$O', 'Mg$_{0.95}$Fe$_{0.05}$O', 'MgO','Mg$_{0.95}$Fe$_{0.05}$O', 'MgO', 'Mg$_{0.98}$Fe$_{0.02}$O','Mg$_{0.95}$Fe$_{0.05}$O', 'Mg$_{0.98}$Fe$_{0.02}$O','Mg$_{0.98}$Fe$_{0.02}$O']
#colors = ['orange','green', 'blue', 'skyblue', 'gold', 'royalblue','coral' , 'navy']


filenames = ['../omegaMgO_39878/39878_Us_T.txt', '../omegaMgO_39874/39874_Us_T.txt', '../omegaMgO_38691/38691_Us_T.txt', '../omegaMgFeO_38692/38692_Us_T.txt', '../omegaMgFeO_39879/39879_Us_T.txt', '../omegaMgFeO_38694/38694_Us_T.txt',    '../omegaMgFeO_38693/38693_Us_T.txt',  '../omegaMgFeO_39877/39877_Us_T.txt',]



colors = ['gold','coral','orange','green', 'limegreen', 'skyblue' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)','Mg$_{0.95}$Fe$_{0.05}$O (38694)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)',  'Mg$_{0.95}$Fe$_{0.05}$O (39877)'  ]




 
   # Soubiran & Militzer 
SMdata_P = [556.425, 631.8, 691.3, 754.38]  #in GPa
SMdata_R = [1.656, 3.017, 4.568, 6.63]      #reflectivity                                                #Burkhardt
#556.425 1.656
#631.8 3.017
#691.3 4.568
#754.38 6.63
    
                                                               
# Reading Lars students file

file_DFT = 'mgo_data_v2.txt'



def combine_files(filename, xx_combined,yy_combined):
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            xx_combined.append(float(columns[0])) #shock velocity
            yy_combined.append(float(columns[3])) #temperature
            
            
    return xx_combined, yy_combined

                                                               
 
def open_files(filename, ax,labelname,color):

    xx = []
    yy = []
    error = []    

    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            xx.append(float(columns[0])) #shock velocity
            yy.append(float(columns[3])) #temperature
        

    x_fit = np.linspace(min(xx),max(xx),100)
    model = np.poly1d(np.polyfit(xx, yy, 2))
    error = np.sqrt(np.sum((model(xx) - yy)**2/len(xx)))


    ax.scatter(xx,yy, s=4, label=labelname,facecolor = color,zorder=2)
    ax.plot(x_fit, model(x_fit),color=color,zorder=2)
    ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2,zorder=2)


x = np.array([15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])

#rho0 = 3.584e3
#a = 1.22568
#b = 7.1552 * 1000

rho0 = 3.584e3
a = 1.2350
b = 7.093172 * 1000

Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments

SMdata_Us = []
f = interpolate.interp1d(Press, x)
for i in range(0,len(SMdata_P)):
    SMdata_Us.append(f(SMdata_P[i]))
    




phase_file = 'phase.txt'
Phase_P = []
Phase_T = []
Phase_Us = []
with open(phase_file, 'r') as file:
    for line in file:
        pressure, temp = line.strip().split(',')
        if 344.0 < float(pressure) :
            Phase_P.append(float(pressure))
            Phase_T.append(float(temp)*1000.0)
            Phase_Us.append(f(pressure))


phase_file = 'phase_B1.txt'
Phase_P_B1 = []
Phase_T_B1 = []
Phase_Us_B1 = []
with open(phase_file, 'r') as file:
    for line in file:
        pressure, temp = line.strip().split(',')
        if 344.0 < float(pressure):
            Phase_P_B1.append(float(pressure))
            Phase_T_B1.append(float(temp)*1000.0)
            Phase_Us_B1.append(f(pressure))

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


V_Mc = []
V_Mc_error = []
T_Mc = []
T_Mc_error = []    

with open('McCoy_R.txt', 'r') as file:
    file.readline()
    for line in file:
        columns = line.strip().split()
        #print(columns)
        V_Mc.append(float(columns[1]))
        V_Mc_error.append(float(columns[2]))
        T_Mc.append(float(columns[5]))
        T_Mc_error.append(float(columns[6]))
        #print(columns[3])
        


T_Mc=np.array(T_Mc)*1000.0
T_Mc_error=np.array(T_Mc_error)*1000.0




fig, ax1 = plt.subplots()


x2 = np.array([16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])

ax1.set_xticks(x2)

#for i in range(0,len(filenames)):
for i in range(0,len(filenames)):
    open_files(filenames[i], ax1,labelnames[i],colors[i])
    
xx_combined=[]
yy_combined=[]
for i in range(0,len(filenames)):
    combine_files(filenames[i],xx_combined,yy_combined)


x_fit = np.linspace(min(xx_combined),max(xx_combined),100)
model = np.poly1d(np.polyfit(xx_combined, yy_combined, 2))
error = np.sqrt(np.sum((model(xx_combined) - yy_combined)**2/len(xx_combined)))

coefficients = model.coefficients
print(coefficients)

#ax1.scatter(xx,yy, s=4, label=labelname,facecolor = color)
ax1.plot(x_fit, model(x_fit),color='grey',alpha=0.2,zorder=1,label='Fit',linewidth='0.5')
ax1.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor='grey', alpha=0.1,zorder=1)

ax1.errorbar(V_Mc, T_Mc, xerr = V_Mc_error, yerr=T_Mc_error, fmt='none', color='tan', label='MgO, M (2019)', elinewidth=0.9)

file_P_T = 'P_T.txt'
f = open(file_P_T, 'w')

xx_combined = np.array(xx_combined)

Press_combined = rho0 * xx_combined * 1000 * (xx_combined * 1000.0-b)/a * 1e-9 #Us-P relationship based on our experiments
for i in range(0,len(xx_combined)):
    f.write(f"{Press_combined[i]} {yy_combined[i]}\n")
f.close()



V_DFT = []
Temp_DFT = []
#Ref_DFT_error = []    

with open(file_DFT, 'r') as file:
    file.readline()
    for line in file:
        columns = line.strip().split()
        V_DFT.append(float(columns[4]))
        Temp_DFT.append(float(columns[2]))
#        Ref_DFT_error.append(float(columns[8]))
        

#errorbar(a, b, yerr=c,
#Ref_DFT=np.array(Ref_DFT)*100.0
#Ref_DFT_error=np.array(Ref_DFT_error)*100.0

ax1.errorbar(V_DFT[0], Temp_DFT[0], fmt='x', color='blue', label='MgO, DFT(B2)')
ax1.errorbar(V_DFT[1:5], Temp_DFT[1:5], fmt='x', color='red', label='MgO, DFT(Liquid)')

#ax1.scatter(SMdata_Us, SMdata_R, label='Soubiran & Militzer (liquid) (2018)')
ax1.plot(Us_smooth, T_smooth,linewidth=0.5,color='steelblue')
ax1.fill_between(Phase_Us,Phase_T, 60000, color='seashell',zorder=0)
ax1.fill_between(Phase_Us,0, Phase_T,color='aliceblue',zorder=0)
ax1.plot(B1_Us_smooth, B1_temperature_smooth,linewidth=0.5,color='gold')
ax1.fill_between(B1_Us_smooth,0, B1_temperature_smooth,color='cornsilk',zorder=0)
ax1.text(15.6, 2000, 'B1 (MgO)',color='gold',fontsize=10)
ax1.text(19, 5000, 'B2 (MgO)',color='steelblue',fontsize=10)
ax1.text(24, 20000, 'Liquid (MgO)',color='lightcoral',fontsize=10)

ax1.set_xlim([15.5, 28])
ylimit = [0,60000]
#ax1.set_ylim([0, 55000]) 
ax1.set_ylim([ylimit[0],ylimit[1] ])
ax1.text(26.5, ylimit[0]+0.1*(ylimit[1]-ylimit[0]),'(b)' ,fontsize=20)



ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1.get_xticks())


#x2 = np.array([16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])
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
ax1.set_ylabel('Temperature (K)', fontsize=17)
ax2.set_xlabel('Pressure (GPa)', fontsize=17)


#ax1.legend(frameon=False,fontsize=8)
fig.savefig('Fig_Us_T.png',dpi=500)
plt.show()
plt.close()