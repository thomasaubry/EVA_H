'''
      	Translated in Python by Biman Chakraborty, April 2022.
      	Originally written in Matlab by Thomas J. Aubry, June 2019.
      	Department of Geography, University of Cambridge
        E-mail: ta460@cam.ac.uk
      	Please cite the corresponding paper if you use this script
'''

#First, import relevant packages and functions to run the script.

import numpy as np
from parameters import ModelParams
from so2injection_8boxes import so2injection_8boxes
from eightboxequations import eightboxequations
from postproc import postproc
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

# ==========================================================================
# 1) ModelParams class from parameters define all model parameter values
# ==========================================================================

modelparam = ModelParams()


# ==========================================================================
# 2) Define volcanic injection parameters
# ==========================================================================

# The function so2injection_8boxes return a list of mass (inmass) and
# time of injection (intime) to be inputed in the model. The main inputs are
# parameters defining the boundaries of the model box and a path to a file
# containing a list of eruptions and their injection source parameters. The
# distribution of SO2 among the 8 model boxes is governed by injection
# latitude and altitude.


#The file volSO2_Carn2016.xlsx contains an inventory of volcanic SO2
#emissions based on Carn et al. 2016, http://dx.doi.org/10.1016/j.jvolgeores.2016.01.002

filename = 'volSO2_Carn2016.xlsx'
print('Calculating volcanic SO2 injections...')
inmass, intime=so2injection_8boxes(modelparam.h1lim,modelparam.h2lim,modelparam.latlim,filename);

# TO RUN THE MODEL FOR A SINGLE ERUPTION INSTEAD, UNCOMMENT THE LINE BELOW.

#inmass, intime=so2injection_8boxes(modelparam.h1lim,modelparam.h2lim,modelparam.latlim,'singleinjection.xlsx');

'''
You can change parameters such as eruption latitude, height and mass in
the file singleinjection.xlsx to see the model sensitivity to eruption
source parameters. If you want more flexibility to look at model
sensitivity, we recommend that you modify the function so2injection_8boxes
to take directly eruption source parameters as inputs instead of an excel
file containing a list of eruption source parameters.
'''

# ==========================================================================
# 3) Define differential equations integration parameters
# ==========================================================================


tspan = [0.5, 455.5-12] # integration bounds in month after jan 1st 1979
# If you use the singleinjection.xlsx file, we recommend that you shorten
# the duration of integration by uncommenting the line below.
# tspan = [0.5,0.5+5*12];


# Set initial conditions ic, i.e. the mass of sulfate in each box, in Tg S
# at first time step. The set of mass entered below correspond to the steady
# state of the model for default parameters values in parameterfile.m. The
# steady state is different from 0 because of background sulfur injections.
ic = np.array([0.0126,0.0468,0.0152,0.0192,0.0359,0.0218,0.0349,0.0417])





# Create a regularly-spaced monthly time vector to interpolate outputs to
# monthly values

tref=np.arange(tspan[0], tspan[-1]+1)

# Create a vector containing time in years, instead of month after Jan 1st
# 1979
time_yr = tref/12+1979  #If you use a different volcanic SO2 inventory, modify
# the start year (1979) to the start year of your inventory.

print('Running the model...')
# ==========================================================================
# 4) Run the model (solve differential equations)
# ==========================================================================

'''
%Run the ode solver. The function eightboxequations contains the
%system of differential equations, to which volcanic SO2 inputs (intime and
%inmass) and model parameters (modelparam and backinj) are passed. The
%outputs of the ode solver are the list of timestep t_ode and the matrix
%containing the mass of sulfate in Tg S in each of the eight boxes, at each
%timestep.
'''
# Note: the tolerance value below have been tested with model timescales as low as 0.1 months.
sol = solve_ivp(eightboxequations, tspan, ic, args=[inmass,intime,modelparam,modelparam.backinj], rtol=1e-4, atol=1e-8) 

# Interpolate outputs to get a monthly timeseries
SO4mass = PchipInterpolator(sol.t,sol.y.T, axis=0)(tref)

print('Doing postprocessing...')
#==========================================================================
#5) Postprocessing
#==========================================================================
#list of wavelengths at which output are requested, in um
wl_req = np.array([0.380, 0.440, 0.525, 0.550, 0.670, 0.870, 1.020])

# %Run the postprocessing function to get essential model outputs:
# %global mean (gm) stratospheric aerosol optical depth (saod), aerosol effective
# %radius (reff), extinction (ext), single scattering albedo (ssa), scattering asymmetry factor
# %(asy). These outputs are dependent on up to 4 variables: time (time_yr),
# %wavelength (wl_req), latitude (lat) and altitude (alt).

gmsaod, saod, reff, ext, ssa, asy, lat, alt=postproc(SO4mass,modelparam,modelparam.mstar,modelparam.R_reff,wl_req);



# Note: the most time-consuming step in this script is the calculation of
# wavelength-dependent properties in the postproc function. If you are only
# interested in outputs at a particular wavelength and need to run the model
# faster, you should request that particular wavelength only when defining
# wl_req above and/or modify the postproc function to request only outputs
# of interest.


print('Making figures...')
# ==========================================================================
# 6) Plot some outputs
# ==========================================================================


# Figure 1: some global mean timeseries
fig1 = plt.figure(figsize=[12,6])

# a) Global mean stratospheric aerosol optical depth at 3 different wavelengths
plt.subplot(3,1,1)
plt.plot(time_yr,gmsaod[:,0])
plt.plot(time_yr,gmsaod[:,3])
plt.plot(time_yr,gmsaod[:,-1])
plt.xlabel('Time (year)')
plt.ylabel('Global mean SAOD')
plt.legend(labels=['380nm','550nm','1020nm'])

# b) Global mean effective radius
plt.subplot(3,1,2)

# calculate global mean mass-weighted effective radius according to the
# scaling from the companion paper
gmreff=modelparam.R_reff*np.sum(SO4mass,axis=1)**(1/3)
#%plot it
plt.plot(time_yr,gmreff)
plt.xlabel('Time (year)')
plt.ylabel('Aerosol effective \n radius ($\mu$ m)')


# c) Global mean net radiative forcing at top of atmosphere: this variable is
# commonly scaled linearly from the global mean SAOD at 550nm. Here we use
# latest estimates the of scaling factor of -24 W/m2 per unit of SAOD, from
# Schmidt et al 2018 (https://doi.org/10.1029/2018JD028776). Note that the
# range of values reported in this particular paper is +/-3W/m2 depending on
# the period used to do the fit, and the range of values found in the
# literature is even larger.
plt.subplot(3,1,3)
plt.plot(time_yr,-24*gmsaod[:,3])
plt.xlabel('Time (year)')
plt.ylabel('Global mean net \n radiative forcing at top \nof atmosphere ($W/m^2$)')
plt.show()

# Figure 2: SAOD at 550nm (log scale)
fig2 = plt.figure(figsize=[8,6])
plt.contour(time_yr,lat,np.log10(saod[:,:,3]).T,colors="black", linestyles="solid")
plt.contourf(time_yr,lat,np.log10(saod[:,:,3]).T,cmap="hot_r")

plt.colorbar(label='SAOD at 550nm (log)', orientation='vertical')
plt.ylabel('Latitude (${}^{o}N$)')
plt.xlabel('Time (year)')
plt.show()



# Figure 3: Extinction near Equator (log scale)
fig3 = plt.figure(figsize=[8,6])
plt.contour(time_yr,alt,np.log10(np.nanmean(ext[:,np.abs(lat)<=5,:,3],axis=0)).T, colors="black", linestyles="solid")
plt.contourf(time_yr,alt,np.log10(np.nanmean(ext[:,np.abs(lat)<=5,:,3],axis=0)).T, cmap="hot_r")
plt.colorbar(label='$5^{o}S-5^{o}N$ mean extinction at 550nm (log)', orientation='vertical')
plt.ylabel('Altitude (km a.s.l)')
plt.xlabel('Time (year)')
plt.show()


# Figure 4: Scattering asymmetry factor at 550nm and 20km a.s.l.
ialt = np.where(alt==20)[0][0]
fig4 = plt.figure(figsize=[8,6])
plt.contour(time_yr,lat,asy[:,:,ialt,3].T, levels=15, colors="black", linestyles="solid")
plt.contourf(time_yr,lat,asy[:,:,ialt,3].T, cmap="hot_r", levels=15)
plt.colorbar(label='Scattering asymmetry factor at 550nm and 20km a.s.l.', orientation='vertical')
plt.ylabel('Latitude (${}^{o}N$)')
plt.xlabel('Time (year)')
plt.show()

#Et voila!
print('Et voila!')