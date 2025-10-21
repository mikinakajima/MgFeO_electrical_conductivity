These files and data are associated with the paper entitled  
"Electrical conductivities of (Mg,Fe)O at extreme pressures and implications for planetary magma oceans" by Nakajima et al.

## **Python scripts**
- `Hugoniot_Up_Us_P_rho.py`: Calculates Hugoniot relationships for (Mg,Fe)O (Figure 2a, 2b)
- `Us_Reflectivity.py`: Calculates the relationship between shock velocity and reflectivity (Figure 3a)
- `Us_Temperature.py`: Calculates the relationship between shock velocity and temperature (Figure 3b)
- `Us_electrical_conductivity.py`: Calculates the relationship between shock velocity and DC conductivity (Figure 3c)
- `Z_values.py`: Calculates the relationship between shock velocity and Z values (Supplementary Figure 2)
- `sigma_model.py`: Calculates the electrical conductivity model (Supplementary Figure 3)

## **Data folders**
- `Hugoniot/`: Stores data for Hugoniot relation plots (Figure 2a, 2b)
- `MgFeO_shock_data/`: Stores data for reflectivity, temperature, DC conductivity, Z values, and conductivity models (Figure 3a–c, Supplementary Figures 2–3)
- `phase/`: Stores phase data files (Figure 3b)

## **Super-Earth modeling**
The super-Earth modeling is based on [Lherm et al. (2024)](https://doi.org/10.1016/j.pepi.2024.107267).  
Relevant codes are available at [https://doi.org/10.17632/3gj9r8rzx9.1](https://doi.org/10.17632/3gj9r8rzx9.1).

## **Note**
- LLNL’s AnalyzeVISAR code was used to analyze the VISAR and SOP data.  
- `IM_MM_MonteCarlo_2019.ipf` written by Dr. Marius Millot was used for impedance matching.
