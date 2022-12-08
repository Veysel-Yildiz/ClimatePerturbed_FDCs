"""
@author: vyildiz
"""
# Import  the modules to be used from Library
import numpy as np
import math 
from scipy import special
import matplotlib.pyplot as plt
import statistics
from func_FDC import *
 
def postplot(num, M, V, L, os_probability, streamflow, av_multiplier, Q_futures , Nsize, low_percentile, case_to_derive): 

 """ 
   This function  plots 4 figures.  
   The ffirst three Figures show the sampling and calculated sattistical paramaters 
   to show if they match each other. 
   The last figure shows a random 3 years of observed stream flow vs derived streamflow
    - num: the size of sampling
    - M: sampled median values 
    - V: sampled coefficient of variation (Cv) values 
    - L: sampled first percentile values
    - os_probability: the exceedance probability of the streamflow records
    - streamflow: input data (observed discharge)
    - av_multiplier: available set of multipliers
    - Nsize:  size of the time series (input)
    - Q_futures: generated future flows
    - low_percentile: the coefficient of low percentile function 
    - case_plot: mean or median case
 """
 

# Figure 1: Fit KOSUGI MODEL to historical data
# Figure 2: derived FDCs


# Derive streamflow statistics
 Q_m, Q_v, Q_low = streamflow_statistics(Q_futures, low_percentile, num, case_to_derive)
 
# Figure 3: plot sampled vs calculated mean/median values
 plt.plot(Q_m, 'ro', label="Derived")
 plt.plot(M, 'b*', label="Sampled")
 plt.legend(loc="upper right")
 plt.grid() 
 plt.xlabel("Futures")
 plt.ylabel("M")
 plt.savefig('PostProcessor_plots' + '/Fig3-M.png')
 plt.clf()

# Figure 4: plot sampled vs calculated Std/CV values
 plt.plot(Q_v, 'ro', label="Derived")
 plt.plot(V, 'b*', label="Sampled")
 plt.legend(loc="upper right")
 plt.grid()
 plt.xlabel("Futures")
 plt.ylabel("V")
 plt.savefig('PostProcessor_plots' + '/Fig4-V.png') 
 plt.clf()


# Figure 5: plot sampled vs calculated low percentile values
 plt.plot(Q_low, 'ro', label="Derived")
 plt.plot(L, 'b*', label="Sampled")
 plt.legend(loc="upper right")
 plt.grid() 
 plt.xlabel("Futures")
 plt.ylabel("Low Percentile [$m^3/s$]")
 plt.savefig('PostProcessor_plots' + '/Fig5-Low.png') 
 plt.clf()


#Figure 6: Random 3 years of observed stream flow vs derived streamflow
 plt.figure(figsize=(11, 6))
 idplot = np.where((av_multiplier[:,1] > 1.75)  & (av_multiplier[:,0] < 0.75) & (0.5 < av_multiplier[:,0]) ) # find the scenario to plot
 idplot = np.asarray(idplot) # converting tuple into int array 
 if np.size(idplot) == 0:
       idplot = np.where(av_multiplier[:,1] >= 1.75)
       idplot = np.asarray(idplot) # converting tuple into int array 
 idplot = np.min(idplot) # get on of the indices if there is more than one
   
 qplot = Q_futures[:,idplot] # select the future 
 qplot = np.reshape(qplot, (len(os_probability),1)) 
 #plt.plot(streamflow[8765:-1],'r')
 #plt.plot(qplot[8765:-1],c='0.35')
 plt.plot(streamflow[8765:-1],'r', label="Observed Streamflow")
 plt.plot(qplot[8765:-1], label="Derived Streamflow",c='0.35')
 plt.legend(loc="upper right")
 plt.xlabel("Time [Days]")
 plt.ylabel("Discharge [$m^3/s$]")
 plt.grid() 
 plt.xlim(0, len(qplot[8765:-1])+10)
 plt.legend(bbox_to_anchor=(1.05, 1))
 plt.tight_layout()
 plt.savefig('PostProcessor_plots' + '/Fig6-ObservedvsDerived_discharge.png') 
 plt.clf()

    

