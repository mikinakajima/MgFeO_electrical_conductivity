import numpy as np
import math
from scipy.constants import pi, G, hbar, m_e, m_u, k, e, N_A
import sys
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


plt.rcParams['font.family'] = 'Helvetica'

file_Us_sigma = 'Us_sigma.txt'


fig, ax1 = plt.subplots()
         
Us, Sigma, PP = np.loadtxt(file_Us_sigma, unpack=True, usecols=(0, 1, 2))        

            
x_fit = np.linspace(min(PP),max(PP),100)
model0, cov_matrix  = np.polyfit(PP, Sigma, 2, cov=True)
model = np.poly1d(model0)
error = np.sqrt(np.sum((model(PP) - Sigma)**2/len(PP)))
coefficients = model.coefficients


errors = np.sqrt(np.diag(cov_matrix))

Us = np.array(Us)

#a_temp = 426.75
#b_temp = -14158.0
#c_temp = 131180.0

rho0, a_v, b_v = np.loadtxt("Up-Us_coefficients.txt")
#b_v = b_v *1000.0
a_temp, b_temp, c_temp = np.loadtxt("T_fit_coefficients.txt")


TT =  a_temp*Us**2.0+b_temp*Us+ c_temp

x_fit2 = np.linspace(min(PP),max(PP),100)


RR=8.3144598
# Define the model function
def model_func(PT, a, b, c):
    P, T = PT
    return a * np.exp(-(b + c * P) /(RR* T))

# Stack P and T together to pass as one argument to curve_fit
PT_data = np.vstack((PP, TT))

# Use curve_fit to find the best parameters
initial_guess = [1, 1, 1]  # Initial guess for a, b, c
params, covariance = curve_fit(model_func, PT_data, Sigma, p0=initial_guess)
errors = np.sqrt(np.diag(covariance))

# Extract the parameters
a, b, c = params
print(f"Best-fit parameters: a = {a}, b = {b}, c = {c}")
print(errors)



# velocity model
#a_v = 1.23011
#b_v = 7.12744

#rho0=3580.0
up=np.linspace(8,15.5) #km
us_fit = a_v * up + b_v #km
Press_fit = rho0 * up*1000.0 * us_fit*1000.0
Press_fit = Press_fit * 1e-9
temp_fit =  a_temp*us_fit**2.0+b_temp*us_fit+ c_temp


sigma_err = np.sqrt(
    (np.exp(-(b + c * Press_fit) / (RR * temp_fit)) * errors[0])**2 +
    ((-a / (RR * temp_fit)) * np.exp(-(b + c * Press_fit) / (RR * temp_fit)) * errors[1])**2 +
    ((-a * Press_fit / (RR * temp_fit)) * np.exp(-(b + c * Press_fit) / (RR * temp_fit)) * errors[2])**2
)


sigma_fit = a * np.exp(-(b + c * Press_fit) /(RR* temp_fit))


for i in ['top', 'bottom','left','right']:
    ax1.spines[i].set_linewidth(1.5)


y = np.array([0, 5e4, 10e4, 15e4, 20e4, 25e4 ])
y_tick_labels= [f'$0$', f'$5$',f'$10$', f'$15$' , f'$20$', f'$25$']

ax1.set_yticklabels(y_tick_labels)
ax1.set_yticks(y)

ax1.scatter(PP,Sigma,label='MgFeO data', s=4)
ax1.plot(Press_fit,sigma_fit,label='Fit', color='grey', alpha=0.5)
ax1.fill_between(Press_fit, sigma_fit - sigma_err, sigma_fit + sigma_err, color='grey', alpha=0.1, edgecolor='none')


ax1.set_xlabel('Pressure (GPa)', fontsize=17)
ax1.set_ylabel('DC Conductivity, $\sigma$ ($10^4$ S/m)', fontsize=17)


ax1.legend(frameon=False,fontsize=8)
fig.savefig('Fig_exp_fit.pdf',dpi=1000)









