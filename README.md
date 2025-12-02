These files and data are associated with the paper entitled  
"Electrical conductivities of (Mg,Fe)O at extreme pressures and implications for planetary magma oceans" by Nakajima et al.

## **Python scripts**
- `Figure2_Hugoniot.py`: Calculates Hugoniot relationships for (Mg,Fe)O (Figure 2a, 2b)
- `Figure3_R_T_Sigma_and_SupplementaryFigure2_Z.py`: Calculates the relationship the reflectivity, temperature, and DC conductivity as functions of shock velocity and pressure (Figure 3a, b, c) and Z values (Supplementary Figure 2) as functions of shock velocity and pressure
- `SupplementaryFigure3_sigma_model.py`: Calculates the electrical conductivity model (Supplementary Figure 3)

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
