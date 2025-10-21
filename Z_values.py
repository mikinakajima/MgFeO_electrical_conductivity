import numpy as np
import math
from scipy.constants import pi, G, hbar, m_e, m_u, k, e, N_A
import sys
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt



plt.rcParams['font.family'] = 'Helvetica'


filenames = ['MgFeO_shock_data/39878_Rfit.txt', 'MgFeO_shock_data/39874_Rfit.txt', 'MgFeO_shock_data/38691_Rfit.txt', 'MgFeO_shock_data/38692_Rfit.txt', 'MgFeO_shock_data/39879_Rfit.txt', 'MgFeO_shock_data/38694_Rfit.txt', 'MgFeO_shock_data/39882_Rfit.txt',   'MgFeO_shock_data/38693_Rfit.txt',  'MgFeO_shock_data/39877_Rfit.txt']
colors = ['gold','coral','orange','green', '#33BEB7', 'skyblue', 'navy' ,'blue', 'royalblue'  ]
labelnames = ['MgO (39878)', 'MgO (39874)','MgO (38691)','Mg$_{0.98}$Fe$_{0.02}$O (38692)', 'Mg$_{0.98}$Fe$_{0.02}$O (39879)','Mg$_{0.95}$Fe$_{0.05}$O (38694)', 'Mg$_{0.95}$Fe$_{0.05}$O (39882)', 'Mg$_{0.95}$Fe$_{0.05}$O (38693)',  'Mg$_{0.95}$Fe$_{0.05}$O (39877)'  ]


file_DFT = 'Hugoniot/mgo_dft_data.txt' #DFT data
file_Us_sigma = 'Us_sigma.txt'

n0 = 1.743 #initial refractive index of MgO from McWilliams 2012
rho0 = 3.58e3  #3.584 
e_0 = 8.854e-12 #Vacuum permittivity
e_charge=1.602e-19 #electron charge, C, A/s
molmass = 40.0 * 1e-3  #kg/mol, MgO, assuming no dissociation

# velocity model
a = 1.23011
b = 7.12744 * 1000

# temperature model
a_temp = 426.75
b_temp = -14158.0
c_temp = 131180.0

f = open(file_Us_sigma, 'w')

def open_files(filename, labelname, color, ax):

        
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
        #print(i,filename)
        Z_ini = 0.0001
        Z = Z_ini # ionization, free parameter  
        r_nb = 0.0
        z_l = 0.0
        z_h = 0.0
        
        
        if RR[i]<=0.0:#(np.abs(n0-1.0)/np.abs(n0+1.0))**2.0:
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
    
        


    ax.scatter(Us,Z_store, s=4, c=color, marker = 'o', alpha=0.3, edgecolors='none') #label=labelname
    
    for i in range(0,len(Us)):
        f.write(f"{Us[i]} {sigma[i]} {pressure[i]}\n")
        

    if len(Us)>0:
        x_fit = np.linspace(min(Us),max(Us),100)
        model = np.poly1d(np.polyfit(Us, Z_store, 2))
        error = np.sqrt(np.sum((model(Us) - Z_store)**2/len(Us)))
        ax.plot(x_fit, model(x_fit),color=color)
        ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2)

 

def calculate_r_nb(rho,Z, temp):


    nb =  1.598 # from Figure S3 McWilliams et al 2012

    ni = (rho*N_A)/molmass #ion density
    ne = Z*ni #electron density 
    v_fermi = (hbar/m_e) * (3 * pi**2 *ne)**(1.0/3.0) #fermi velocity 
    v_thermal = np.sqrt((2.0*k*temp)/m_e) #thermal velocity
    vel = [v_fermi, v_thermal] # use the larger velocity between fermi and thermal velocities
    l =  2.0 * (3.0/(4.0*np.pi*ni))**(1.0/3.0)

    omega = 2.0 * np.pi * 2.99792e8/(532e-9) #frequency at 532 nm
    omega_p = np.sqrt((ne * e_charge**2)/(e_0*m_e)) #plasma frequency
    scatter_tau = l/max(vel)#v_fermi
    
    omegatau=scatter_tau * omega
    sigma0 = (ne * e_charge **2 * scatter_tau)/m_e #DC conductivity eq. 43 from Millot 2015    
    n_mod =np.sqrt(nb**2 - (omega_p**2) / (omega**2 * (1 + 1j/(omega*scatter_tau))))  
    r_nb = (np.abs(n0-n_mod)/np.abs(n0+n_mod))**2.0
    return r_nb, sigma0, omegatau


fig, ax1 = plt.subplots()
for i in range(0,len(filenames)):
    open_files(filenames[i], labelnames[i],colors[i],ax1)
    
x = np.arange(16,29)
ax1.set_xticks(x)

y = np.array([0, 0.1, 0.2,0.3, 0.4])
ax1.set_yticks(y)

f.close()
ax1.set_xlim([15.5, 28])

ylimit = [0,0.35]
ax1.set_ylim([ylimit[0],ylimit[1]])
Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments
Temp =  a_temp * x**2.0 + b_temp*x + c_temp


# Load DFT data (skip header)
data = np.loadtxt(file_DFT, skiprows=1)
V_DFT = data[:, 4]
Sigma_DFT = data[:, 5]
Sigma_DFT_error = data[:, 6]


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

ax2.set_xticks(velocity_for_press_ticks)
ax2.set_xticklabels(custom_press_ticks)


for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)


ax1.tick_params(axis='both', which='major', labelsize=11)
ax2.tick_params(axis='both', which='major', labelsize=11)


ax1.set_xlabel('Shock Velocity (km/s)', fontsize=16)
ax1.set_ylabel('Z values', fontsize=16)
ax2.set_xlabel('Pressure (GPa)', fontsize=16)

for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)


ax1.legend(frameon=False,fontsize=8, loc = "upper left", bbox_to_anchor=(0.01, 0.95))
fig.savefig('Fig_Z.pdf',dpi=1000)



plt.show()
plt.close()

