import scipy as sci 
from scipy.integrate import solve_ivp
import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import scipy.interpolate
import matplotlib.pyplot as plt
import statsmodels.api as sm
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, ScalarFormatter, MaxNLocator)

'''Loading data from excel'''

df_GWB = pd.read_excel('Fe_Os_GWB_Compilation_Submission.xlsx')
Fe_ultra = df_GWB['f-Ultramafic 1:4']
Ca_SO4_relation_ultra = df_GWB['Ca/SO4']
Fe_bas = df_GWB['f-Basalt 1:4']
Ca_SO4_relation_bas = df_GWB['Ca/SO4']



df_Ca_SO4 = pd.read_excel("Ca_SO4_data_compilation_Submission.xlsx")
age_Ca_SO4 = df_Ca_SO4['Age (Ma)']*1E6
Ca_SO4_record = df_Ca_SO4['Ca/SO4 median']

#################################################################
'''Setting up variables to be used later on in the code. These values can be manipulated to change, for example, the values of the Os mass balance'''

#time step in yrs
t_interval = 1000

#number of resampling
num_monte = 100000

#time length of model to run
t_end = 7E7

#Array holding total run time, split up be time step
t_change1 = np.arange(0, t_end, 1000000)


#Os mass balance parameters, all fluxes in mol/yr

Os_ht = 2.8 #high Temp hydrothermal Os flux
R_Os_ht = 0.1296 #high Temp hydrothermal Os ratio

Os_bw = 580 #island basalt weathering Os flux
R_Os_bw = 0.1296 #island basalt weathering Os ratio

Os_cosm = 80 #cosmogenic Os flux
R_Os_cosm = 0.128 #cosmogenic Os ratio

Os_dust = 36.8 #dust Os flux
R_Os_dust = 1.05 #dust Os ratio

Os_riv = 1800 #Riverine Os flux
R_Os_riv = 1.44 #Riverine Os ratio

Os_lowt = 62 #low Temp hydrothermal Os flux 
R_Os_lowt = 0.87 #low Temp hydrothermal Os ratio

#Amount of Os in ocean in mol
Os_ocean = 7.20e7

#Modern seawater Os ratio
R_SW = 1.06

#Hydrothermal fluid flux (kg/yr), from Elderfield and Shultz, 1996 
HT_flux = 6E13

#Fraction of hydrothermal fluid flux that alters ultramafic- and basalt-hosted systems
frac_ultra = 0.25
frac_bas = 0.75



#################################################################
'''Manipulations to data for smoothing and interpolating, as well as unit conversions'''

#Smoothing and interpolation for SO4_Ca ratios
Ca_SO4_record_smooth = lowess(Ca_SO4_record, age_Ca_SO4, frac = 0.9, missing = "drop")
Ca_SO4_record_interp = sci.interpolate.interp1d(Ca_SO4_record_smooth[:,0], Ca_SO4_record_smooth[:,1], kind = "cubic")

sc_min = lowess(df_Ca_SO4['Ca/SO4 min'], df_Ca_SO4['Age (Ma)'], frac = 0.85)
sc_max = lowess(df_Ca_SO4['Ca/SO4 max'],df_Ca_SO4['Age (Ma)'], frac = 0.85)

sc_min_interp = sci.interpolate.interp1d(sc_min[:,0], sc_min[:,1], kind = 'cubic')
sc_max_interp = sci.interpolate.interp1d(sc_max[:,0], sc_max[:,1], kind = 'cubic')

sc_result_max = sc_max_interp(df_Ca_SO4['Age (Ma)'])
sc_result_min = sc_min_interp(df_Ca_SO4['Age (Ma)'])


#Set up array to hold seawater chemistry data, including error 
nrow = df_Ca_SO4.shape[0] - 1
df_Ca_SO4_smooth = np.zeros(shape=(nrow,num_monte))
for i in range(nrow):
    df_Ca_SO4_smooth[i,:] = np.random.uniform(low=sc_result_min[i], high=sc_result_max[i], size=num_monte)


#Generate array of size num_monte holding random uniform distribution of initial Os (fmol/kg)
Os_ini_monte_ultra = np.random.uniform(1500, 1500000, num_monte)
Os_ini_monte_bas = np.random.uniform(1500, 1500, num_monte)

#Generate array of size num_monte holding random uniform distribution of partition coefficient
dist_coef_monte = np.random.uniform(10, 20, num_monte) 


#Generate empty array to hold results below, lines 112-118
Os_ultra_remain_monte = np.zeros((len(Fe_ultra), num_monte))
Os_bas_remain_monte = np.zeros((len(Fe_bas), num_monte))

#Calculating how much Os left with sulfide formation/Fe change in ultramafic systems, in fmol/kg
for i in range(0, num_monte):
    Os_ultra_remain_monte[:,i] = (((Os_ini_monte_ultra[i]*4)+(50*1))/(4+1)) * Fe_ultra**(dist_coef_monte[i])
    
#Calculating how much Os left with sulfide formation/Fe change in basalt systems, in fmol/kg
for i in range(0, num_monte):
    Os_bas_remain_monte[:,i] = (((Os_ini_monte_bas[i]*4)+(50*1))/(4+1)) *  Fe_bas**(dist_coef_monte[i])


#flux of Os remaining in fluid after sulfide formation/Fe change in mol/yr
Os_ultra_remain_flux = (Os_ultra_remain_monte *1E-15) * HT_flux * frac_ultra
Os_bas_remain_flux = (Os_bas_remain_monte *1E-15) * HT_flux * frac_bas

#Set up empty lists to collect data from below
Os_ini_suc_ultra = []
Os_ini_suc_bas = []
dist_coef_suc = []
HT_flux_bal_record = []
Os_remain_flux_suc = pd.DataFrame()
dOs_plot_monte = pd.DataFrame()


t = 0
z = 0

for i in range(0, num_monte):
    Os_Ca_SO4_relation_interp_ultra = sci.interpolate.interp1d(Ca_SO4_relation_ultra, Os_ultra_remain_flux[:,i])
    
    Os_Ca_SO4_relation_interp_bas = sci.interpolate.interp1d(Ca_SO4_relation_bas, Os_bas_remain_flux[:,i])
        
    #Setting the initial values for various variables
    dOs_ocean_old = np.array(1.06/(7.4 + 1.06))
    
    dOs_ocean = dOs_ocean_old
    
    
    #Use Sharma et al., 2007 imbalance calculation to filter for mass balance reasons
    HT_flux_bal = Os_Ca_SO4_relation_interp_ultra(Ca_SO4_record_interp(t_change1[0])) + Os_Ca_SO4_relation_interp_bas(Ca_SO4_record_interp(t_change1[0]))
    HT_flux_bal_record.append(HT_flux_bal)
    
    if (HT_flux_bal * (R_Os_ht - 1.06) + Os_dust * (R_Os_dust - 1.06) + Os_cosm * (R_Os_cosm - 1.06) + Os_bw * (R_Os_bw - 1.06) + Os_riv * (R_Os_riv - 1.06) + Os_lowt * (R_Os_lowt - 1.06)) >= -1 and (HT_flux_bal * (R_Os_ht - 1.06) + Os_dust * (R_Os_dust - 1.06) + Os_cosm * (R_Os_cosm - 1.06) + Os_bw * (R_Os_bw - 1.06) + Os_riv * (R_Os_riv - 1.06) + Os_lowt * (R_Os_lowt - 1.06)) <= 1:
        
        Os_ini_suc_ultra.append(Os_ini_monte_ultra[i])
        Os_ini_suc_bas.append(Os_ini_monte_bas[i])
        dist_coef_suc.append(dist_coef_monte[i])
        Os_remain_flux_suc.insert(len(Os_remain_flux_suc.columns), z, Os_ultra_remain_flux[:,i] + Os_bas_remain_flux[:,i])
        z = z + 1


#normalizing ratios to get more accurate calculation (see Li and Elderfield 2013)
R_Os_ht = R_Os_ht / (7.4 + R_Os_ht)
R_Os_lowt = R_Os_lowt / (7.4 + R_Os_lowt)
R_Os_cosm = R_Os_cosm / (7.4 + R_Os_cosm)
R_Os_dust = R_Os_dust / (7.4 + R_Os_dust)
R_Os_bw = R_Os_bw / (7.4 + R_Os_bw)
R_Os_riv = R_Os_riv / (7.4 + R_Os_riv)
R_Os_modern_SW = R_SW / (7.4 + R_SW)

#Set up and run the ODE to solve the forward Os mass balance model through time
z = 0 
for i in range(0, len(Os_remain_flux_suc.columns)):
    lowess = sm.nonparametric.lowess
    w = lowess(df_Ca_SO4_smooth[:,i], df_Ca_SO4['Age (Ma)'][0:6]*1E6, frac=0.9)
    
    Ca_SO4_fun = sci.interpolate.interp1d(w[:,0], w[:,1])
    
    Os_Ca_SO4_relation_interp = sci.interpolate.interp1d(Ca_SO4_relation_bas, Os_remain_flux_suc[i])
    
    K_Os_flux_times_ratio = lambda t, y: Os_Ca_SO4_relation_interp(Ca_SO4_fun(t)) * (R_Os_ht - y) + Os_bw * (R_Os_bw - y) + Os_cosm * (R_Os_cosm - y) + Os_dust * (R_Os_dust - y) + Os_riv * (R_Os_riv - y) + Os_lowt * (R_Os_lowt - y)
        
    K_Os_slope_ratio = lambda t, y: K_Os_flux_times_ratio(t,y)/Os_ocean
        
    K_Os_ode_result = solve_ivp(K_Os_slope_ratio, (0, t_end), dOs_ocean_old.flatten(), method = 'LSODA', t_eval = t_change1)
    
    dOs_total = K_Os_ode_result.y
    
    t_plot = t_change1
    
    dOs_plot = dOs_total.flatten()
    dOs_plot = (dOs_plot * 7.4)/(1 - dOs_plot)
    dOs_plot_monte.insert(len(dOs_plot_monte.columns), z ,dOs_plot)
    
    z = z+1

#Calculate mean and standard deviation of the monte carlo routine
dOs_plot_monte_mean = np.mean(dOs_plot_monte,1)
dOs_plot_monte_std = np.std(dOs_plot_monte,1)


#Age manipulation to make plots
age_model = t_plot/1E6
#################################################################

'''Plotting results'''
f, ax3 = plt.subplots(figsize=(7.09, 7.09))
ax3.plot(age_model, dOs_plot_monte_mean, color = ('#00356B'), label = 'Model')
ax3.fill_between(age_model, dOs_plot_monte_mean + dOs_plot_monte_std, dOs_plot_monte_mean - dOs_plot_monte_std, facecolor= ('#00356B'), alpha = 0.2)
ax3.set_ylabel('$^{187}$Os/$^{188}$Os',fontsize = 16, fontname = 'Arial',  labelpad = 10)
ax3.set_xlabel('Age (Ma)',fontsize = 16, fontname = 'Arial',  labelpad = 10)
ax3.tick_params(axis = 'both', length = 5, width = 1.5, labelsize = 14, direction = 'in', pad = 10)
ax3.tick_params(which = 'minor', length = 3, width = 1.2, direction = 'in')
ax3.set_xlim(-0.5,68)
ax3.set_ylim(0.1,1.3)
ax3.xaxis.set_major_locator(MultipleLocator(10))
ax3.xaxis.set_major_formatter(ScalarFormatter())
ax3.xaxis.set_minor_locator(MultipleLocator(5))
