#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:51:28 2023

@author: mikinakajima
"""

# This python script plots Figure 2 of Nakajima et al.
# Impedance matching is not conducted in this script. 
# Instead, it was conducted using IM_MM_MonteCarlo_2019.ipf written by Dr. Marius Millot.


import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import sys
from scipy.optimize import curve_fit

plt.rcParams['font.family'] = 'Helvetica'



rho0=3.580 #MgO density


colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue', 'navy' ,'blue', 'royalblue']
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','(Mg$_{0.98}$,Fe$_{0.02}$)O (38692)', '(Mg$_{0.98}$,Fe$_{0.02}$)O (39879)','(Mg$_{0.95}$,Fe$_{0.05}$)O (38694)', '(Mg$_{0.95}$,Fe$_{0.05}$)O (39882)', '(Mg$_{0.95}$,Fe$_{0.05}$)O (38693)',  '(Mg$_{0.95}$,Fe$_{0.05}$)O (39877)'  ]
labelnames_composition = ['MgO', 'MgO','MgO','(Mg$_{0.98}$,Fe$_{0.02}$)O', '(Mg$_{0.98}$,Fe$_{0.02}$)O ','(Mg$_{0.95}$,Fe$_{0.05}$)O', '(Mg$_{0.95}$,Fe$_{0.05}$)O ', '(Mg$_{0.95}$,Fe$_{0.05}$)O',  '(Mg$_{0.95}$,Fe$_{0.05}$)O'  ]
labelnames_shotID = ['39878', '39874','38691','38692', '39879','38694', '39882', '38693',  '39877'  ]
thickness = [76, 76, 62, 58, 76, 58, 77, 58, 80]
Energy = [529.4, 595.7, 784.5, 809.4, 613.8, 365.2, 640.0, 731.0, 828.4]



def linear_model(x, a, b):
    return a * x + b

def poly_model(x, a, b, c):
    return a * x**2.0 + b * x + c



def load_block(path, cols, skiprows=1):
    """Generic loader: path + column map (0-based indices)."""
    use = [c for c in cols.values() if c is not None]
    arr = np.loadtxt(path, skiprows=skiprows, usecols=use, unpack=False)
    if arr.ndim == 1:
        arr = arr[None, :]
    pos = {c: i for i, c in enumerate(use)}
    return {k: (arr[:, pos[c]] if c is not None else None) for k, c in cols.items()}

# column maps for each dataset
MCCOY       = {"us":3,"up":5,"us_err":4,"up_err":6,"p":7,"p_err":8,"rho":9,"rho_err":10}
ROOT        = {"us":3,"up":1,"us_err":4,"up_err":2,"p":7,"p_err":8,"rho":5,"rho_err":6}
MCWILLIAMS  = {"us":3,"up":5,"us_err":4,"up_err":6,"p":7,"p_err":8,"rho":9,"rho_err":10}
OURS        = {"us":6,"up":8,"us_err":7,"up_err":9,"p":2,"p_err":3,"rho":4,"rho_err":5}

# DFT-MD calculations: rho=1, P=3, V=4 (0-based indices)
DFT = {"rho":1, "p":3, "V":4}

# load datasets
mccoy = load_block("Hugoniot/McCoy.txt", MCCOY)
root  = load_block("Hugoniot/Root.txt", ROOT)
mcw   = load_block("Hugoniot/McWilliams.txt", MCWILLIAMS)
ours  = load_block("Hugoniot/MgO_ours.txt", OURS)
dft   = load_block("Hugoniot/MgO_dft_data.txt", DFT)   
        
# combining all the previous + our data    
us_all = np.concatenate([mccoy["us"], root["us"], mcw["us"], ours["us"]])
up_all = np.concatenate([mccoy["up"], root["up"], mcw["up"], ours["up"]])
p_all = np.concatenate([mccoy["p"], root["p"], mcw["p"], ours["p"]])
rho_all = np.concatenate([mccoy["rho"], root["rho"], mcw["rho"], ours["rho"]])


# Up-Us relationship fit
params, covariance = curve_fit(linear_model, up_all, us_all)
a, b = params


a_err = np.sqrt(covariance[0, 0])  # Standard error of 'a'
b_err = np.sqrt(covariance[1, 1])  # Standard error of 'b'
print(a,a_err,b, b_err, 'a,b, Us-Up')

y_param = np.linspace(5.0, 20,100)
x_param = a*y_param + b

mse = np.mean((linear_model(up_all,a,b) - us_all) ** 2)


plt.figure(figsize=(12, 15))
ax = plt.subplot(211)  

plt.xticks(fontsize=25)
plt.yticks(fontsize=25)


plt.yticks(np.arange(10, 35, 5))

plt.errorbar(mccoy["up"], mccoy["us"], xerr=mccoy["up_err"], yerr=mccoy["us_err"], fmt='o', markersize = 3, color='brown',label ='MgO (McCoy et al. 2019)')
plt.errorbar(root["up"], root["us"], xerr=root["up_err"], yerr=root["us_err"], fmt='o', markersize = 3, color='grey', label ='MgO (Root et al. 2019)')
plt.errorbar(mcw["up"], mcw["us"], xerr=mcw["up_err"], yerr=mcw["us_err"], fmt='o', markersize = 3, color='purple', label ='MgO (McWilliams et al. 2012)')



for i in ['top','bottom','left','right']:
    ax.spines[i].set_linewidth(3)



for i in range(0,9):
    plt.errorbar(ours["up"][i], ours["us"][i], xerr=ours["up_err"][i], yerr=ours["us_err"][i], fmt='o', markersize = 15, color = colors[i], label=labelnames[i])#, label ='This work')

plt.plot(y_param, x_param, color='black', label='Fit')


plt.xlabel('Particle Velocity (km/s)', fontsize=30)
plt.ylabel('Shock Velocity (km/s)', fontsize=30)
plt.legend(frameon=False,fontsize=13,loc='upper left')

plt.text(0.9, 0.05, "(a)", transform=plt.gca().transAxes, fontsize=25, fontweight='bold')


ax2 = plt.subplot(212)  # 1 row, 2 columns, first subplot

plt.text(0.9, 0.05, "(b)", transform=plt.gca().transAxes, fontsize=25, fontweight='bold')

for i in ['top','bottom','left','right']:
    ax2.spines[i].set_linewidth(3)

plt.errorbar(mccoy["rho"], mccoy["p"], xerr=mccoy["rho_err"], yerr=mccoy["p_err"], fmt='o', color='brown', markersize = 3)#, label ='MgO (McCoy et al. 2019)')
plt.errorbar(root["rho"], root["p"], xerr=root["rho_err"], yerr=root["p_err"],color='grey', fmt='o', markersize = 3)#, label ='MgO (Root et al. 2019)')
plt.errorbar(mcw["rho"], mcw["p"], xerr=mcw["rho_err"], yerr=mcw["p_err"], color='purple',fmt='o', markersize = 3)#, label ='MgO (McWilliams et al. 2012)')

for i in range(0,9):
  plt.errorbar(ours["rho"][i], ours["p"][i], xerr=ours["rho_err"][i], yerr=ours["p_err"][i], fmt='o', markersize = 15, color = colors[i])# label ='MgO, MgFeO (this work)')
 

plt.errorbar(dft["rho"][0], dft["p"][0], fmt='x', markersize = 10, color='blue', label ='MgO,DFT-MD(B2)')
plt.errorbar(dft["rho"][1:5], dft["p"][1:5], fmt='x', markersize = 10, color='red', label ='MgO,DFT-MD(B2)')



#rho-P relationship fit
params2, covariance = curve_fit(poly_model, rho_all, p_all)
a2, b2, c2 = params2

a2_err = np.sqrt(covariance[0, 0])  
b2_err = np.sqrt(covariance[1, 1])  
c2_err = np.sqrt(covariance[2, 2])  

print(a2,a2_err,b2, b2_err, c2, c2_err, 'a,b,c, rho-P')


x_param2 = np.linspace(5, 10,100)
y_param2 = a2*x_param2**2.0 + b2* x_param2 + c2
plt.plot(x_param2, y_param2, color='black', label='Fit')


x_axis = [5, 6, 7, 8, 9, 10] 
plt.xticks(x_axis)


for i in range(len(ours["us"])):
    print(
        f"{labelnames_shotID[i]} & {labelnames_composition[i]} & {thickness[i]} & {Energy[i]:.1f} "
        f"& {ours['p'][i]:.2f} $\\pm$ {int(ours['p_err'][i])} "
        f"& {ours['rho'][i]:.3f} $\\pm$ {ours['rho_err'][i]:.3f} "
        f"& {ours['us'][i]:.2f} $\\pm$ {ours['us_err'][i]:.2f} "
        f"& {ours['up'][i]:.2f} $\\pm$ {ours['up_err'][i]:.2f} \\\\"
    )


plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim([5.5, 10])
plt.legend(frameon=False,fontsize=13,loc='upper left')

plt.xlabel('Density (g cm$^{-1}$)', fontsize=30)
plt.ylabel('Pressure (GPa)', fontsize=30)

plt.savefig('Fig_Us_Up_rho_P.eps',dpi=1000)
plt.show()
plt.close()
