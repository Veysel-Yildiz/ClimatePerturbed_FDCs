"""
@author: veysel yildiz
"""

# Import  the modules to be used from Library
import numpy as np
import pandas as pd
import math 
from scipy import special
import scipy.optimize
import matplotlib.pyplot as plt
from smt.sampling_methods import LHS
from operator import itemgetter
import statistics
from scipy.optimize import root

# Import  the all the functions defined
from func_FDC import *
from PostProcessor_FDC import *


""" 
FIRST SECTION: Input Section ..................................................1
    
"""

# Load the input data set
streamflow = np.loadtxt('input' + '/b_observed' + '.txt', dtype=float, delimiter=',')


###Select the methodological case; the key streamflow characteristics
# Median Case: Median, Coefficient of Variation and Low Percentile (first or fifth)
# Mean Case: Mean, standard deviation and Low Percentile (first or fifth)
    
methodological_case = 'Median' # Type the case name : Median / Mean
low_percentile = 1 # Define low flow quantile: 1 or 5

# sample size of future scenarios
num = 50

#Define scaling factos ranges for pars:  M, V, L
xlimits = np.array([[0.3, 1], [1.0, 2.0], [0.3, 1]]) 

######################################################################################
""" 
SECOND SECTION: Fit a FDC to historical records................................2
"""

## 1.1 Derive  FDC from historical data
fdc_discharge, fdc_probability  = Calc_Fdc(streamflow);  # Derive FDC

## 1.2 Fit Kosugi Model to historical data
FDC_pars, p_pred, K_Rmse = Fdc_Metrics(fdc_discharge,fdc_probability); # Calculate the FDC metrics for the fitted curve

# Derive discharge from Kosugi Model
k_discharge = Kosugi(FDC_pars, p_pred)

##
if methodological_case == 'Mean':
 case_to_derive = 1
  # define the functions for STD
 Der_Analytical = Std_analytical # analytical solution for std
 Der_Opt = Std_Opt # define pars with optimisation  for std
 

if methodological_case == 'Median':
 case_to_derive = 2 
 # define the functions for STD
 Der_Analytical = CV_analytical # analytical solution for std
 Der_Opt = CV_Opt # define pars with optimisation  for std
 
E = math.exp(math.sqrt(2)*special.erfcinv(2*(1- low_percentile/100))) # calculate the coefficient of low percentile function    
 
 #Convert a 1D array to a 2D Numpy array for streamflow_statistics function 
streamflow_2d = streamflow.reshape(len(streamflow), 1)
##

 # Derive streamflow statistics  of historical records
M1, V1, L1 = streamflow_statistics(streamflow_2d, low_percentile, 1, case_to_derive)
 
#####################################################################################

#derive parameter FDC parameters of historical records
FDC_pars_m1  = Der_Analytical(M1,V1,L1,E); 

# Derive discharge for  multiplier of 1
m1_discharge = Kosugi(FDC_pars_m1, fdc_probability)

#Calculate RMSE 
#M1_Rmse = round(Opt_Rmse(FDC_pars_m1, fdc_discharge,fdc_probability)[0],4)

## 1.4: Test the applicability of the approach

# Figure 1: plot historical FDC vs best fit and curve of parameters of historical records
plt.rcParams["figure.figsize"] = (15,8)
plt.plot(fdc_probability*100, fdc_discharge,marker='o', markeredgecolor='r', linestyle='none', markerfacecolor='none', label="Historical Data")
plt.plot(p_pred*100, k_discharge,'k', label="Best Fit")
plt.plot(fdc_probability*100, m1_discharge,'b', label="FDC from $M_h, V_h, L_h$")
#plt.plot(p_pred*100, k_discharge,'k', label="Fitted Curve, RMSE $="  + str(K_Rmse) + "$")
#plt.plot(fdc_probability*100, m1_discharge,'b', label="Multiplier of 1, RMSE $="  + str(M1_Rmse) + "$")
plt.legend(loc="upper right")
plt.xlabel("Exceedance Probability [%]", labelpad=20)
plt.ylabel("Flow rate [$m^3/s$]", labelpad=20)
plt.rc('font', size=20)
plt.grid() 
plt.savefig('PostProcessor_plots' + '/Fig1-Fitted FDC.png')
plt.clf()

# Clear unnecessary variables from memory
del   K_Rmse, p_pred,  FDC_pars, FDC_pars_m1, streamflow_2d 

"""
THIRD SECTION: LHS Sampling for statistical parameters.........................3
"""

sampling = LHS(xlimits=xlimits)

# Use ranges and sampling to get multiplier for all three inputs
multiplier = sampling(num) 

M = multiplier[:, 0]*M1  #Mean Sampling

V = multiplier[:, 1]*V1# Variance Sampling

L = multiplier[:, 2]*L1 #First percntile Sampling

# find the scenarios in which sampled mean is bigger than  sampled first percentile
if L1 > 0:
 diff = M / L 
 idx_s = np.where(diff > 1.3)  
else : 
 diff = M - L
 idx_s = np.where(diff > 0.3) 

idx = np.transpose(np.asarray(idx_s)) # get the index number of the available scenarios
av_num = len(idx); # find the number of available scenarios

#select the available scenarios
M = M[idx]
V = V[idx]
L = L[idx]

#available set of multipliers
new_multip= multiplier[idx, :]
av_multiplier = np.reshape(new_multip, (av_num,3))

# Clear unnecessary variables from memory
del num, idx_s, idx, M1, L1, V1, new_multip, diff, multiplier, xlimits

"""
FOURTH SECTION: Generate ensemble of future scenarios..........................4
"""

Nsize = np.size(fdc_probability) # find the size of the time series (input)
Q_futures = np.empty((Nsize,av_num)) # for future flows
Q_futures [:] = np.NaN
b = np.empty((av_num,1)) 
b [:] = np.NaN

#return exceedance probability of the originial sequence of the streamflow
os_probability = Return_OS(streamflow, fdc_probability);  
 
    
counter=0 # adding a counter
# Iterates over the scenarios to derive new FDC pars: a_s, b_s, c_s
for i in range(av_num):
  
 der_par = lambda par: Der_Opt( M[i], V[i], L[i], E, fdc_probability, Nsize, par)[0] #define the function to be optimized
 b[i] = scipy.optimize.fmin(func=der_par, x0=2, disp=False) # do the optimization
 
 Fdc_pars_future = Der_Opt( M[i], V[i], L[i], E, fdc_probability, Nsize, b[i])[1]
 Q_futures[:,i] = Kosugi(Fdc_pars_future, os_probability)
 Q_futures[:,i][Q_futures[:,i] < 0.01] = 0   
 counter+=1
 print(counter)    
    

# Figure: derived FDCs
for j in range(av_num):
    s_fdc = Calc_Fdc(Q_futures [:,j]);  # Derive FDC
    s_probability = s_fdc[1] #Define the exceedance probability
    s_discharge = s_fdc[0] #Define the streamflow
    plt.plot(s_probability*100, s_discharge, c='0.35')
    plt.yscale("log")
    if j == av_num-1:
        l_fdc = Calc_Fdc(Q_futures[:,av_num-1]);  # Derive FDC
        plt.plot(l_fdc[1]*100, l_fdc[0], label="Derived FDCs",c='0.35')
        plt.plot(fdc_probability*100, fdc_discharge,'b', label="Streamflow Records")
        plt.legend(loc="upper right")
        plt.xlabel("Exceedance Probability [%]", labelpad=20)
        plt.ylabel("Log Flow rate [$m^3/s$]", labelpad=20)
        plt.rcParams["figure.figsize"] = (10,5)
        plt.rc('font', size=20)
        plt.grid() 
        plt.ylim(0.01, 500)
        plt.savefig('PostProcessor_plots' + '/Fig2-Future FDCs.png')
        plt.clf()


# Clear unnecessary variables from memory
del   i, j, counter,  s_fdc, s_discharge, s_probability,  l_fdc, Fdc_pars_future, b

"""
FIFTH SECTION: Post processor .................................................5
"""
Post_Processor = 'Yes' # Type 'Yes' to plot the figures or Type 'No' 

if Post_Processor == 'Yes':
 postplot(av_num, M, V, L, os_probability, streamflow, av_multiplier, Q_futures, Nsize, low_percentile, case_to_derive)

# Clear unnecessary variables from memory
#del av_num, Post_Processor, av_multiplier, Nsize, E

