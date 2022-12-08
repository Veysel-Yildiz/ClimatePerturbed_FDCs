"""@author: vyildiz
"""
# Import  the modules to be used from Library
import scipy.optimize
import numpy as np
import math 
from scipy import special
#from scipy.optimize import fsolve
from scipy.optimize import root
from operator import itemgetter
import statistics

##################################################11111111111111111111111111111
def Calc_Fdc(streamflow):

 """ 
    This function is used to derive a FDC which is  a plot of discharge vs. percent of time 
    that a particular discharge was equaled or exceeded.
    - streamflow: observed streamflow
 """
 
# sorting in descending order
 x = np.sort(streamflow)[::-1]

# Create a a column vector  for probability
 nl = np.arange(1.,len(x)+1)

# calculate the probabilities of the streamflow
 p = (nl-0.5)/ len(x) 

# Now return the exceedance probability
 return x,p

##################################################111111111111111111111111111111
def Kosugi(FDC_pars, p_pred):

 """ 
    Kosugi model to derive  streamflow
    - a,b,c: fitting parameters of Kosugi Model
    - p_pred: exceedance probability
 """
 
 # Unpack  pars:a,b,c of Kosugi Model
 a = FDC_pars[0]
 b = FDC_pars[1]
 c = FDC_pars[2]
 
 #p_pred[0] = 1.75e-05
 def f(p_pred):
    return c + (a-c)*np.exp(math.sqrt(2)*b*special.erfcinv(2*p_pred))
 f_vectorize = np.vectorize(f)
 k_discharge = f_vectorize(p_pred)

# Now return the discharge
 return k_discharge

##################################################222222222222222222222222222222
def Return_OS(streamflow,fdc_probability):

 """ 
 This function returns exceedance probability of the originial sequence of the streamflow
    - streamflow: observed streamflow
    - fdc_probability:  the exceedance probability of historical records
 """
 
# sort the streamflow ank keep  the list of indices for the  originial sequence
 indices, q_sorted = zip(*sorted(enumerate(-streamflow), key=itemgetter(1)))

 idx_sorted = np.asarray(indices)
 
#return the originial sequence of the streamflow
 idx = np.empty(len(fdc_probability)) #  create a new array for originial sequence

 for i in range(len(fdc_probability)): 
  idn = np.where(idx_sorted==i)
  idx[i] = np.asarray(idn)

 idx = idx.astype(int) # converting the float array into int array 

# make sure that it returns to the original sequence
 #streamflow_sorted  = -np.asarray(q_sorted)
 #original_sequence_streamflow = fdc_discharge[idx] 

 f_pred = fdc_probability[idx] 
 
# exceedance probability of  originial sequence
 return f_pred 

##################################################333333333333333333333333333333
def Fdc_Metrics(x,p): 

 """ 
    This function is used to calculate the Fdc metrics of the discharge data
    - x : the exceedance probability
    - p : the streamflow
 """
 
 # define the function to be optimized
 f_pars = lambda pars: Opt_Rmse(pars,x,p)[0]
 X_opt = scipy.optimize.fmin(func=f_pars, x0=[0.5,1.5,0.1],disp=False) # do the optimization

 p_pred = Opt_Rmse(X_opt,x,p)[1]
 Rmse = round(Opt_Rmse(X_opt,x,p)[0],4)
 
 return X_opt, p_pred, Rmse

#####################################################444444444444444444444444444
def Opt_Rmse(pars,x,p): 

 """ 
    This function is used to calculate RMSE of the proposed parameters
    - pars: FDC pars; a,b and c
    - x : the exceedance probability
    - p : the streamflow
 """

# Unpack pars
 a = pars[0]
 b = pars[1]
 c = pars[2]
 
 #create an array with a specific size for probability
 p_pred = np.zeros(len(x)) 
 
 #predict the probabilities
 for i in range(len(x)):
     if x[i]> c:
       p_pred[i] = math.erfc(1/(math.sqrt(2)*b)*math.log((x[i]-c)/(a-c)))/2
     else:
      p_pred[i] = 1
      
 # Calculate RMSE
 RMSE = ( sum(( p_pred - p)**2) )**0.5 / ( len(x))

 return RMSE, p_pred

###########################################################555555555555555555555

def CV_analytical(M, V, L, E): 

 """ 
    Root Function to derive FDC parameter-b values 
    
    - M: sampled median values 
    - V: sampled coefficient of variation (Cv) values 
    - L: sampled low percentile values
    - E: the coefficient of low percentile function 
 """


 solve_CV = root(lambda b: math.sqrt(np.exp(b**2) - 1) / (((L / M - E**b)/(1-L / M)) /  np.exp(b**2/2) + 1)-V, V)
 
 b = solve_CV.x;
 
  # Derive new pars based on par b
 a = M

 c = (L- M * E**b) / (1 - E**b )
 
 return np.array([a, b, c])

###########################################################666666666666666666666


def Std_analytical(M, V, L, E): 

 """ 
    Root Function to derive FDC parameter-b values 
    
    - M: sampled mean values 
    - V: sampled standard deviation values
    - P: sampled low percentile values
    - E: the coefficient of low percentile function 
 """

 solve_Var = root(lambda b: ((M - L) / (np.exp(b**2/2) - E**b)) * (np.exp(2*b**2) - np.exp(b**2))**0.5-V, 2)
 b = solve_Var.x;
 
 # Derive new pars based on par b
 a = ( L*(np.exp(b**2/2) - 1) + M* (1- E**b) ) / ( (np.exp(b**2/2) - E**b) )
 
 c = ( L* np.exp(b**2/2) - M * E**b)/ ( (np.exp(b**2/2) - E**b) )
 
 return np.array([a, b, c])

###########################################################777777777777777777777
def CV_Opt(M, V, L, E,e,N_size, b): 

 """ 
    Root Function to derive FDC parameter-b values 
    
    - M: sampled median values 
    - V: sampled coefficient of variation (Cv) values 
    - P: sampled first percentile values
    - E: the coefficient of low percentile function 
 """

   # Define f(b)
  #vectorized version of the function: to accept a single element to every element in an array
 def f(e):
    return np.exp(math.sqrt(2)*b*special.erfcinv(2*e))
 f1 = np.vectorize(f)
 f_b = sum(f1(e))/N_size # f(b)

  # Define f(2b)
 def f(e):
    return (np.exp(math.sqrt(2)*b*special.erfcinv(2*e)))**2
 f2= np.vectorize(f)
 f_2b = sum(f2(e))/N_size # f(2b)

 a = M
 c = (L- M * E**b) / (1 - E**b ) 
 
 Q_v = (a -c) * math.sqrt(f_2b - f_b**2) / (c + (a-c)* f_b)

#Error function to be minimized
 Rerr =  (Q_v - V)**2
 fdcpars = np.array([a, b, c])
 
 return Rerr, fdcpars

###########################################################888888888888888888888

def Std_Opt(M, V, L, E,e,N_size, b): 

 """ 
    Root Function to derive FDC parameter-b values 
    
    - M: sampled mean values 
    - V: sampled standard deviation values
    - P: sampled low percentile values
    - E: the coefficient of low percentile function 
 """

   # Define f(b)
  #vectorized version of the function: to accept a single element to every element in an array
 def f(e):
    return np.exp(math.sqrt(2)*b*special.erfcinv(2*e))
 f1 = np.vectorize(f)
 f_b = sum(f1(e))/N_size # f(b)

  # Define f(2b)
 def f(e):
    return (np.exp(math.sqrt(2)*b*special.erfcinv(2*e)))**2
 f2= np.vectorize(f)
 f_2b = sum(f2(e))/N_size # f(2b)

 c = ( L* f_b - M * E**b)/  (f_b - E**b) 
 a = ( L*(f_b - 1) + M* (1- E**b) ) / ( f_b - E**b) 
 
 Q_v = (a-c)*math.sqrt(f_2b - f_b**2)

#Error function to be minimized
 Rerr =  (Q_v - V)**2
 fdcpars = np.array([a, b, c])
 
 return Rerr, fdcpars


 
###########################################################888888888888888888888

def streamflow_statistics(Q_futures, low_percentile, num, case_to_derive): 

 """ 
    Derive streamflow statistics of future FDCs
    
    - Q_futures: derived future FDCs
    - low_percentile: the coefficient of low percentile function 
    - num: number of futures
    - case_to_derive: mean or median case
 """
 
 Q_m = np.empty((num)) #  create a new array for median
 Q_v = np.empty((num)) #  create a new array for CV
 Q_low = np.empty((num)) #  create a new array for first percentile
 #Q_m[:] = np.NaN
   
 if case_to_derive == 1 : #'Mean'
    for i in range(num):  
        Q_v[i] = Q_futures[:,i].std()  # calculate V
        Q_m[i] = statistics.mean(Q_futures[:,i]) # calculate mean
        Q_low[i]  = np.percentile(Q_futures[:,i], low_percentile) #calculate the first percntile
     
 elif case_to_derive == 2:
     for i in range(num):  
         Q_v[i] = Q_futures[:,i].std() /  Q_futures[:,i].mean()# calculate V
         Q_m[i] = statistics.median(Q_futures[:,i]) # calculate mean
         Q_low[i]  = np.percentile(Q_futures[:,i], low_percentile) #calculate the first percntile
     
 return  Q_m, Q_v, Q_low


