#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:51:28 2023

@author: mikinakajima
"""

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from scipy.optimize import curve_fit

plt.rcParams['font.family'] = 'Helvetica'

x_values = []
y_values = []
x_error = []
y_error = []
Press = []
Press_error = []
rho = []
rho_error = []
rho0=3.580

yvalue = 5 # 7: pressure, 5: Particle velocity I guess?

with open('McCoy.txt', 'r') as file:
    next(file)
    for line in file:

        columns = line.strip().split()
        x_values.append(float(columns[3])) #shock velocity
        y_values.append(float(columns[yvalue]))
        x_error.append(float(columns[4]))
        y_error.append(float(columns[yvalue+1]))

        Press.append(float(columns[7]))
        Press_error.append(float(columns[8]))        
        rho.append(float(columns[9]))
        rho_error.append(float(columns[10]))       
        
        
        
        
x_root = []
y_root = []
x_root_error = []
y_root_error = []
Press_root = []
Press_root_error = []
rho_root = []
rho_root_error = []



yvalue = 1

with open('Root.txt', 'r') as file:
    next(file)
    for line in file:
        columns = line.strip().split()
        x_root.append(float(columns[3]))
        y_root.append(float(columns[yvalue]))
        x_root_error.append(float(columns[4]))
        y_root_error.append(float(columns[yvalue+1]))
        
       
        Press_root.append(float(columns[7]))
        Press_root_error.append(float(columns[8]))        
        rho_root.append(float(columns[5]))
        rho_root_error.append(float(columns[6]))  
       
        
      
x = []
y = []
x_e = []
y_e = []

Press_MW = []
Press_MW_error = []
rho_MW = []
rho_MW_error = []


yvalue = 5

with open('McWilliams.txt', 'r') as file:
    next(file)
    for line in file:

        columns = line.strip().split()
        x.append(float(columns[3]))
        y.append(float(columns[yvalue]))
        x_e.append(float(columns[4]))
        y_e.append(float(columns[yvalue+1]))
        
        Press_MW.append(float(columns[7]))
        Press_MW_error.append(float(columns[8]))        
        rho_MW.append(float(columns[9]))
        rho_MW_error.append(float(columns[10]))  
        
        
x2 = []
y2 = []
x2_e = []
y2_e = []
Press2 = []
Press2_error = []
rho2 = []
rho2_error = []
rho_mod =[]
Press_mod=[]
rho_error_mod =[]
Press_error_mod=[]

# old data
#39882 '(Mg$_{0.95}$,Fe$_{0.05}$)O' 1030.51 31.53 8.625 0.461 22.18 0.42 12.96 0.374 0.570
#39878 'MgO' 1143.92 30.39 8.708 0.386 23.29 0.31 13.70 0.348 0.574

#38692 '(Mg$_{0.98}$,Fe$_{0.02}$)O'  1261.91 37.31 8.75 0.43 24.42 0.33 14.41 0.42 0.
#38694 '(Mg$_{0.95}$,Fe$_{0.05}$)O'  761.61 29.20  8.05 0.44 19.57 0.38 10.85 0.39 0.532 

yvalue = 8

#Alex - you can just change MgO.txt to add your data points.

with open('MgO_v2.txt', 'r') as file:
    next(file)
    for line in file:

        columns = line.strip().split()
        x2.append(float(columns[6])) #shock velocity
        y2.append(float(columns[yvalue])) #particle velocity
        x2_e.append(float(columns[7])) #shock velocity error
        y2_e.append(float(columns[yvalue+1])) #particle velocity error
        
        Press2.append(float(columns[2]))
        Press2_error.append(float(columns[3]))    
        
        Us=float(columns[6])
        Up=float(columns[8])
        dUs=float(columns[7])
        dUp=float(columns[9])
        
        rho_mod.append(rho0* Us/(Us-Up))
        Press_mod.append(rho0* Us*Up)  
        
        # Compute dP using error propagation
        dP = rho0* Us*Up * np.sqrt((dUs / Us) ** 2 + (dUp / Up) ** 2)

# Compute dRho using error propagation
        dRho = np.sqrt((rho0 * dUs / (Us - Up)) ** 2 + (rho0 * Us * dUp / (Us - Up) ** 2) ** 2)

        rho_error_mod.append(dRho)
        Press_error_mod.append(dP)  
        
        rho2.append(float(columns[4]))
        rho2_error.append(float(columns[5]))  
        
        
        
x_fit = np.concatenate((x_values, x_root, x, x2))
y_fit = np.concatenate((y_values, y_root, y, y2))
def linear_model(x, a, b):
    return a * x + b


V_DFT=[]
P_DFT=[]
rho_DFT=[]


with open('mgo_data_revised.txt', 'r') as file:
    next(file)
    for line in file:

        columns = line.strip().split()
        V_DFT.append(float(columns[4]))
        P_DFT.append(float(columns[3]))
        rho_DFT.append(float(columns[1]))



params, covariance = curve_fit(linear_model, y_fit, x_fit)
a, b = params

a_error = np.sqrt(covariance[0, 0])  # Standard error of 'a'
b_error = np.sqrt(covariance[1, 1])  # Standard error of 'b'

y_param = np.linspace(5.0, 20,100)
x_param = a*y_param + b

mse = np.mean((linear_model(y_fit,a,b) - x_fit) ** 2)


plt.figure(figsize=(12, 15))
ax = plt.subplot(211)  # 1 row, 2 columns, first subplot

plt.xticks(fontsize=25)
plt.yticks(fontsize=25)


plt.yticks(np.arange(10, 35, 5))

plt.errorbar(y_values, x_values, xerr=y_error, yerr=x_error, fmt='o', markersize = 3, color='brown',label ='MgO (McCoy et al. 2019)')
plt.errorbar(y_root, x_root, xerr=y_root_error, yerr=x_root_error, fmt='o', markersize = 3, color='grey', label ='MgO (Root et al. 2019)')
plt.errorbar(y, x, xerr=y_e, yerr=x_e, fmt='o', markersize = 3, color='purple', label ='MgO (McWilliams et al. 2012)')



for i in ['top','bottom','left','right']:
    ax.spines[i].set_linewidth(3)

#colors = ['orange','green','blue', 'skyblue','gold','royalblue', 'coral', 'limegreen', 'navy' ]
#labelnames = ['MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)', 'Mg$_{0.95}$Fe$_{0.05}$O (38694)','MgO (39874)', 'Mg$_{0.95}$Fe$_{0.05}$O (39877)', 'MgO (39878)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)',  'Mg$_{0.95}$Fe$_{0.05}$O (39882)']


colors = ['gold','coral','orange','green', 'limegreen', 'skyblue', 'navy' ,'blue', 'royalblue',  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)','Mg$_{0.95}$Fe$_{0.05}$O (38694)', 'Mg$_{0.95}$Fe$_{0.05}$O (39882)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)',  'Mg$_{0.95}$Fe$_{0.05}$O (39877)'  ]

labelnames_composition = ['MgO', 'MgO','MgO','Mg$_{0.98}$Fe$_{0.02}$O', 'Mg$_{0.98}$Fe$_{0.02}$O ','Mg$_{0.95}$Fe$_{0.05}$O', 'Mg$_{0.95}$Fe$_{0.05}$O ', 'Mg$_{0.95}$Fe$_{0.05}$O',  'Mg$_{0.95}$Fe$_{0.05}$O'  ]
labelnames_shotID = ['39878', '39874','38691','38692', '39879','38694', '39882', '38693',  '39877'  ]
thickness = [76, 76, 62, 58, 76, 58, 77, 58, 80]
Energy = [529.4, 595.7, 784.5, 809.4, 613.8, 365.2, 640.0, 731.0, 828.4]


for i in range(0,9):
    plt.errorbar(y2[i], x2[i], xerr=y2_e[i], yerr=x2_e[i], fmt='o', markersize = 15, color = colors[i], label=labelnames[i])#, label ='This work')



plt.plot(y_param, x_param, color='black', label='Fit')
print(a,b,a_error,b_error, 'a,b, param')

plt.xlabel('Particle Velocity (km/s)', fontsize=30)
plt.ylabel('Shock Velocity (km/s)', fontsize=30)
plt.legend(frameon=False,fontsize=13,loc='upper left')

plt.text(0.9, 0.05, "(a)", transform=plt.gca().transAxes, fontsize=25, fontweight='bold')


ax2 = plt.subplot(212)  # 1 row, 2 columns, first subplot

plt.text(0.9, 0.05, "(b)", transform=plt.gca().transAxes, fontsize=25, fontweight='bold')

for i in ['top','bottom','left','right']:
    ax2.spines[i].set_linewidth(3)

plt.errorbar(rho, Press, xerr=rho_error, yerr=Press_error, fmt='o', color='brown', markersize = 3)#, label ='MgO (McCoy et al. 2019)')
plt.errorbar(rho_root, Press_root, xerr=rho_root_error, yerr=Press_root_error,color='grey', fmt='o', markersize = 3)#, label ='MgO (Root et al. 2019)')
plt.errorbar(rho_MW, Press_MW, xerr=rho_MW_error, yerr=Press_MW_error, color='purple',fmt='o', markersize = 3)#, label ='MgO (McWilliams et al. 2012)')

for i in range(0,9):
 #   plt.errorbar(rho2[i], Press2[i], xerr=rho2_error[i], yerr=Press2_error[i], fmt='o', markersize = 15, color = colors[i])# label ='MgO, MgFeO (this work)')
     plt.errorbar(rho_mod[i], Press_mod[i], xerr=rho_error_mod[i], yerr=Press_error_mod[i], fmt='o', markersize = 15, color = colors[i])# label ='MgO, MgFeO (this work)')



#plt.errorbar(x_values, Press, xerr=x_error, yerr=Press_error, fmt='o', color='brown', markersize = 3)# label ='MgO (McCoy et al. 2019)')
#plt.errorbar(x_root, Press_root, xerr=x_root_error, yerr=Press_root_error,color='grey', fmt='o', markersize = 3)# label ='MgO (Root et al. 2019)')
#plt.errorbar(x, Press_MW, xerr=x_e, yerr=Press_MW_error, color='purple',fmt='o', markersize = 3)#, label ='MgO (McWilliams et al. 2012)')


plt.errorbar(rho_DFT[0], P_DFT[0], fmt='x', markersize = 10, color='blue', label ='MgO,DFT-MD(B2)')
plt.errorbar(rho_DFT[1:5], P_DFT[1:5], fmt='x', markersize = 10, color='red', label ='MgO,DFT-MD(liquid)')


#for i in range(0,9):    
#    plt.errorbar(x2[i], Press2[i], xerr=x2_e[i], yerr=Press2_error[i], fmt='o', markersize = 15, color = colors[i])# label ='MgO, MgFeO (this work)')

x_fit2 = np.concatenate((rho, rho_root, rho_MW, rho2))
y_fit2 = np.concatenate((Press, Press_root, Press_MW, Press2))

def poly_model(x, a, b, c):
    return a * x**2.0 + b * x + c

params2, covariance = curve_fit(poly_model, x_fit2, y_fit2)
a2, b2, c2 = params2

a_2 = np.sqrt(covariance[0, 0])  # Standard error of 'a'
b_2 = np.sqrt(covariance[1, 1])  # Standard error of 'b'
c_2 = np.sqrt(covariance[2, 2])  # Standard error of 'c'
#b_2 = np.sqrt(covariance[1, 1])  # Standard error of 'b'



x_param2 = np.linspace(5, 10,100)
y_param2 = a2*x_param2**2.0 + b2* x_param2 + c2
plt.plot(x_param2, y_param2, color='black', label='Fit')


x_axis = [5, 6, 7, 8, 9, 10]  # Generates 7 points between 16 and 28
plt.xticks(x_axis)

for i in range(len(x2)):
    print(f"{labelnames_shotID[i]} & {labelnames_composition[i]} & {thickness[i]} & {Energy[i]:.1f} "
          f"& {Press_mod[i]:.2f} $\\pm$ {int(Press_error_mod[i])} "
          f"& {rho_mod[i]:.3f} $\\pm$ {rho_error_mod[i]:.3f} "
          f"& {x2[i]:.2f} $\\pm$ {x2_e[i]:.2f} "
          f"& {y2[i]:.2f} $\\pm$ {y2_e[i]:.2f} \\\\")


plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim([5.5, 10])
plt.legend(frameon=False,fontsize=13,loc='upper left')

plt.xlabel('Density (g cm$^{-1}$)', fontsize=30)
plt.ylabel('Pressure (GPa)', fontsize=30)

plt.savefig('Fig_Us_Up_rho_P.eps',dpi=500)
plt.show()
plt.close()