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
rho0 = 3.58e3  
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

    Us = []
    RR = []
    error = []  
    rho = []
    up = []
    temp = []
    pressure = []

    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            
            Us.append(float(columns[1]))
            RR.append(float(columns[0])*0.01)


    for i in range(0,len(Us)):
        up.append((Us[i] * 1000 -b)/a)
        rho.append(rho0 * Us[i]*1000/(Us[i]*1000-up[i])) # density after shock (at Us= 25km/s)
        temp.append((a_temp * Us[i]**2.0 + b_temp*Us[i] + c_temp))
        pressure.append((rho0 * Us[i]*1000* up[i])*1e-9)
        

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
        #print(omegatau0)
        


    ax.scatter(Us,sigma, s=4, c=color, marker = 'o', alpha=0.3, edgecolors='none') #label=labelname
    
    for i in range(0,len(Us)):
        f.write(f"{Us[i]} {sigma[i]} {pressure[i]}\n")
        

    if len(Us)>0:
        x_fit = np.linspace(min(Us),max(Us),100)
        model = np.poly1d(np.polyfit(Us, sigma, 2))
        error = np.sqrt(np.sum((model(Us) - sigma)**2/len(Us)))
        ax.plot(x_fit, model(x_fit),color=color)
        ax.fill_between(x_fit, model(x_fit)+ error,  model(x_fit)- error, facecolor=color, alpha=0.2)

 

def calculate_r_nb(rho,Z, temp):


    nb =  1.598 # from Figure S3 McWilliams et al 2012

    ni = (rho*N_A)/molmass #ion density
    ne = Z*ni #electron density 

    #temp = 30000
    v_fermi = (hbar/m_e) * (3 * pi**2 *ne)**(1.0/3.0) #fermi velocity 
    v_thermal = np.sqrt((2.0*k*temp)/m_e) #thermal velocity
    vel = [v_fermi, v_thermal] # use the larger velocity between fermi and thermal velocities
    l =  2.0 * (3.0/(4.0*np.pi*ni))**(1.0/3.0)

    omega = 2.0 * np.pi * 2.99792e8/(532e-9) 
    #omega = 2.99e8/(532e-9) #frequency at 532 nm
    omega_p = np.sqrt((ne * e_charge**2)/(e_0*m_e)) #plasma frequency
    scatter_tau = l/max(vel)#v_fermi
    
    omegatau=scatter_tau * omega
    
    
    #scatter_tau_org = (2.0/max(vel)) * (3.0/(4.0*pi**2.0 *ni))**(0.333) #tau scatter time, s #pi seems wrong!
    sigma0 = (ne * e_charge **2 * scatter_tau)/m_e #DC conductivity eq. 43 from Millot 2015

    n_mod =np.sqrt(nb**2.0-omega_p**2.0/omega**2.0/(1+1j/omega/scatter_tau))
    r_nb = (np.abs(n0-n_mod)/np.abs(n0+n_mod))**2.0

    #n_real = np.sqrt(nb**2.0 - (omega_p/omega)**2.0/(1.0 + complex(0,1) /omega/scatter_tau)).real #real part of index of refraction
    #n_imag = np.sqrt(nb**2.0 - (omega_p/omega)**2.0/(1.0 + complex(0,1) /omega/scatter_tau)).imag #real part of index of refraction
    #r_nb = ((n_real - n0)**2.0 +n_imag**2.0)/ ((n_real + n0)**2 + n_imag**2.0) #reflectivity 

    return r_nb, sigma0, omegatau
#    print('reflectivity = ', r_nb, 'DC conductivity (S/m)=',sigma0, 'Carrier Density (1/cm^3)=',ne*1e-6)


fig, ax1 = plt.subplots()
for i in range(0,len(filenames)):
    open_files(filenames[i], labelnames[i],colors[i],ax1)
    

x = np.array([16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])
ax1.set_xticks(x)
#ax1.set_yscale("log")

y = np.array([0, 5e4, 10e4, 15e4, 20e4, 25e4, 30e4 ])
y_tick_labels= [f'$0$', f'$5$',f'$10$', f'$15$' , f'$20$', f'$25$', f'$30$']

ax1.set_yticklabels(y_tick_labels)
ax1.set_yticks(y)


f.close()

#ax1.set_ylim([10, 250000])
ax1.set_xlim([15.5, 28])


ylimit = [-10000,250000]
ax1.set_ylim([ylimit[0],ylimit[1]])
ax1.text(26.5, ylimit[0]+0.1*(ylimit[1]-ylimit[0]),'(c)' ,fontsize=20)

Press = rho0 * x * 1000 * (x * 1000-b)/a * 1e-9 #Us-P relationship based on our experiments
Temp =  a_temp * x**2.0 + b_temp*x + c_temp



V_DFT = []
Sigma_DFT = []
Sigma_DFT_error = []    

with open(file_DFT, 'r') as file:
    file.readline()
    for line in file:
        columns = line.strip().split()
        V_DFT.append(float(columns[4]))
        Sigma_DFT.append(float(columns[5]))
        Sigma_DFT_error.append(float(columns[6]))

Sigma_DFT=np.array(Sigma_DFT)
Sigma_DFT_error=np.array(Sigma_DFT_error)

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
from scipy.interpolate import interp1d
pressure_to_velocity = interp1d(Press, x, bounds_error=False, fill_value="extrapolate")
velocity_for_press_ticks = pressure_to_velocity(custom_press_ticks)

# Set ticks and labels on the top axis
ax2.set_xticks(velocity_for_press_ticks)
ax2.set_xticklabels(custom_press_ticks)





# strixrude comparison
sigma00 = 1.99e5

Delta_E = 75.94e3
Delta_V = -0.061 * 1e-6
#P = 500e9
R = 8.3
TT = 40000.0 #this needs to be updated I think

P = Press * 1e9
sigma_stix=np.zeros(len(x))
for i in range(0,len(x)):
    sigma_stix[i] = sigma00 * np.exp(-(Delta_E+P[i]*Delta_V)/(R*Temp[i]))

ax1.plot(x, sigma_stix*10**(-4), label='MgO ', color = 'tan', linestyle='-')
ax1.plot(x, sigma_stix*10**(-3.2), label='Mg$_{0.95}$Fe$_{0.05}$O ', color = 'navy', linestyle='dashed')
ax1.plot(x, sigma_stix, label='Mg$_{0.75}$Fe$_{0.25}$O ', color = 'grey', linestyle='dashed')
ax1.text(15.7, 23.5*1e4, 'DFT-MD (Holmström et al. 2018)', fontsize=8)

ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)

#def combine_files(filename, xx_combined,yy_combined):
xx_combined=[]
yy_combined=[]
with open(file_Us_sigma, 'r') as file:
    for line in file:
        columns = line.strip().split()
        xx_combined.append(float(columns[0])) #shock velocity
        yy_combined.append(float(columns[1])) #reflectivity
            
        #print(columns[0], columns[1])
    
xx_mod=[]#xx_combined
yy_mod=[]#yy_combined

xx_mod_s=[]
yy_mod_s=[]

for i in range(0,len(xx_combined)):
    if xx_combined[i]>0: #19.5
        xx_mod.append(xx_combined[i])
        yy_mod.append(yy_combined[i])  
    else:
        xx_mod_s.append(xx_combined[i])
        yy_mod_s.append(yy_combined[i])  
    



x_fit = np.linspace(min(xx_mod),max(xx_mod),100)
model0, cov_matrix = np.polyfit(xx_mod, yy_mod, 2, cov=True) #np.poly1d(np.polyfit(xx_mod, yy_mod, 2 ))
model = np.poly1d(model0)
error = np.sqrt(np.sum((model(xx_mod) - yy_mod)**2/len(xx_mod)))

coefficients = model.coefficients
errors = np.sqrt(np.diag(cov_matrix))


x_fit2 = np.linspace(min(xx_mod),max(xx_mod),100)
#ax1.plot(x_fit2,coefficients[0]*x_fit2**2.0 + coefficients[1]*x_fit2 +coefficients[2] ,color='grey',alpha=1,zorder=0,linewidth='1')
#ax1.fill_between(x_fit2, model(x_fit2)+ error,  model(x_fit2)- error, facecolor='grey', alpha=0.1,zorder=0)

print("Coefficients:", coefficients)
print("errors:", errors)


for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)






ax1.set_xlabel('Shock Velocity (km/s)', fontsize=17)

ax1.set_ylabel('DC Conductivity, $\sigma$ ($10^4$S/m)', fontsize=17)


ax1.legend(frameon=False,fontsize=8, loc = "upper left", bbox_to_anchor=(0.01, 0.95))
fig.savefig('Fig_sigma.pdf',dpi=1000,  bbox_inches="tight")



plt.show()
plt.close()

