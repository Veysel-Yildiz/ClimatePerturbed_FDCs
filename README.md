
This repository contains Python code, along with a Jupyter Notebook demo, for perturbing the flow duration curve to explore future flow conditions.
The methodology presented in the paper by  V. Yildiz, C. Rougé, R. Milton, S.F. Brown "Technical note: Statistical generation of climate-perturbed flow duration curves" submitted to the Hydrology and Earth System Sciences.

The library is built exclusively using Python code; It consists of three files. 

Contents:

`Der_FDC.py`: This is the main file to run. It consists of five sections. 
Section-1 requires users to place the historical streamflow records in the input file and to select the case to derive futures; "median" or "mean". What is more, users can change number of scenarios and the range for statistical parameters in this part.  

Section-2 derives flow duration curve (FDC) from the historical data provided. Also it fits Kosugi Model to historical data and curve of multiplier of 1 of the parameters. The code  generates a plot of historical FDC vs fitted curves on (Figure-1).  

Section-3 samples future scenarios by using the concept of latin hypercube sampling (LHS). 
Section-4 derives new FDC pars (a,b andc) values based on sampled scenarios and then generate future flows 
The second figure is plotted shows the the future FDCs based on sampling. 
Finally, in Section-5, post processor figures are plotted and saved in PostProcessor_plots file.  

`func_FDC.py`: This file contains all the  functions which are used in the  main file (Der_FDC). It is advised not to modify this file. 

`PostProcessor_FDC.py`: This function plots 4 figures. 
Fig.3, 4,5 show  the sampling and calculated sattistical paramaters(mean/median, standard deviation/coefficient of variation, low percentile) to show if they match each other or not. The last figure (Fig. 6) shows a 3 year subset derived streamflows based on a random future versus 3 years of observed streamflow values for the same time period (last three years)..

