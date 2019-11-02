% 	Written by Thomas J. Aubry, June 2019.
% 	Department of Geography, University of Cambridge
%   E-mail: ta460@cam.ac.uk
% 	Please cite the corresponding paper if you use this script

clear
close all

%==========================================================================
%1) Run parameterfile.m to define all model parameter values
%==========================================================================
parameterfile


%==========================================================================
%2) Define volcanic injection parameters
%==========================================================================

%The function so2injection_8boxes return a list of mass (inmass) and
%time of injection (intime) to be inputed in the model. The main inputs are
%parameters defining the boundaries of the model box and a path to a file
%containing a list of eruptions and their injection source parameters. The
%distribution of SO2 among the 8 model boxes is governed by injection
%latitude and altitude.


%The file volSO2_Carn2016.xlsx contains an inventory of volcanic SO2
%emissions based on Carn et al. 2016, http://dx.doi.org/10.1016/j.jvolgeores.2016.01.002
disp('Calculating volcanic SO2 injections...')
[inmass intime]=so2injection_8boxes(h1lim,h2lim,latlim,'volSO2_Carn2016.xlsx');

%TO RUN THE MODEL FOR A SINGLE ERUPTION INSTEAD, UNCOMMENT THE LINE BELOW.

 %[inmass intime]=so2injection_8boxes(h1lim,h2lim,latlim,'singleinjection.xlsx');

%You can change parameters such as eruption latitude, height and mass in
%the file singleinjection.xlsx to see the model sensitivity to eruption
%source parameters. If you want more flexibility to look at model
%sensitivity, we recommend that you modify the function so2injection_8boxes
%to take directly eruption source parameters as inputs instead of an excel
%file containing a list of eruption source parameters.

%==========================================================================
%3) Define differential equations integration parameters
%==========================================================================


tspan = [0.5 455.5-12];%integration bounds in month after jan 1st 1979
%If you use the singleinjection.xlsx file, we recommend that you shorten
%the duration of integration by uncommenting the line below.
% tspan = [0.5 0.5+5*12];


%Set initial conditions ic, i.e. the mass of sulfate in each box, in Tg S
%at first time step. The set of mass entered below correspond to the steady
%state of the model for default parameters values in parameterfile.m. The
%steady state is different from 0 because of background sulfur injections.
ic = [0.0126;0.0468;0.0152;0.0192;0.0359;0.0218;0.0349;0.0417];


%Set tolerance for the differential equation solver. The value below have
%been tested with model timescales as low as 0.1 months.
ode_opts = odeset('RelTol',1e-4,'AbsTol',1e-8);


%Create a regularly-spaced monthly time vector to interpolate outputs to
%monthly values
tref=(tspan(1):1:tspan(end))';

%Create a vector containing time in years, instead of month after Jan 1st
%1979
time_yr=tref/12+1979;%If you use a different volcanic SO2 inventory, modify
%the start year (1979) to the start year of your inventory.

disp('Running the model...')
%==========================================================================
%4) Run the model (solve differential equations)
%==========================================================================

%Run the ode45 solver in matlab. See
%https://uk.mathworks.com/help/matlab/ref/ode45.html for more details on
%this Matlab's ode solver. The function eightboxequations contains the
%system of differential equations, to which volcanic SO2 inputs (intime and
%inmass) and model parameters (modelpara and backinj) are passed. The
%outputs of the ode solver are the list of timestep t_ode and the matrix
%containing the mass of sulfate in Tg S in each of the eight boxes, at each
%timestep.
[t_ode,SO4mass] = ode45(@(t,y) eightboxequations(t,y,inmass,intime,modelpara,backinj), tspan, ic, ode_opts);%

%Interpolate outputs to get a monthly timeseries
SO4mass=interp1(t_ode,SO4mass,tref,'pchip');

disp('Doing postprocessing...')
%==========================================================================
%5) Postprocessing
%==========================================================================
%list of wavelengths at which output are requested, in um
wl_req=[0.380 0.440 0.525 0.550 0.670 0.870 1.020]';

%Run the postprocessing function to get essential model outputs:
%global mean (gm) stratospheric aerosol optical depth (saod), aerosol effective
%radius (reff), extinction (ext), single scattering albedo (ssa), scattering asymmetry factor
%(asy). These outputs are dependent on up to 4 variables: time (time_yr),
%wavelength (wl_req), latitude (lat) and altitude (alt).

[gmsaod, saod, reff, ext, ssa, asy, lat, alt]=postproc(SO4mass,modelpara,mstar,R_reff,wl_req);



%Note: the most time-consuming step in this script is the calculation of
%wavelength-dependent properties in the postproc function. If you are only
%interested in outputs at a particular wavelength and need to run the model
%faster, you should request that particular wavelength only when defining
%wl_req above and/or modify the postproc function to request only outputs
%of interest.


disp('Making figures...')
%==========================================================================
%6) Plot some outputs
%==========================================================================


%Figure 1: some global mean timeseries
figure(1)

%a) Global mean stratospheric aerosol optical depth at 3 different wavelengths
subplot(3,1,1)
plot(time_yr,gmsaod(:,1))
hold on
plot(time_yr,gmsaod(:,4))
plot(time_yr,gmsaod(:,end))
xlabel('Time (year)')
ylabel('Global mean SAOD')
legend('380nm','550nm','1020nm')

%b) Global mean effective radius
subplot(3,1,2)

%calculate global mean mass-weighted effective radius according to the
%scaling from the companion paper
gmreff=R_reff*sum(SO4mass,2).^(1/3);
%plot it
plot(time_yr,gmreff)
xlabel('Time (year)')
ylabel({'Aerosol effective';'radius (\mu m)'})


%c) Global mean net radiative forcing at top of atmosphere: this variable is
%commonly scaled linearly from the global mean SAOD at 550nm. Here we use
%latest estimates the of scaling factor of -24 W/m2 per unit of SAOD, from
%Schmidt et al 2018 (https://doi.org/10.1029/2018JD028776). Note that the
%range of values reported in this particular paper is +/-3W/m2 depending on
%the period used to do the fit, and the range of values found in the
%literature is even larger.
subplot(3,1,3)
plot(time_yr,-24*gmsaod(:,4))
xlabel('Time (year)')
ylabel({'Global mean net';'radiative forcing at top';'of atmosphere (W/m2)'})


%Figure 2: SAOD at 550nm (log scale)
figure(2)
contourf(time_yr,lat,squeeze(log10(saod(:,:,4)))')
colormap(flipud(hot))
h=colorbar;
ylabel(h,'SAOD at 550nm (log)')
ylabel('Latitude (^{o}N)')
xlabel('Time (year)')



%Figure 3: Extinction near Equator (log scale)
figure(3)
contourf(time_yr,alt,squeeze(log10(nanmean(ext(:,abs(lat)<=5,:,4),2)))')
colormap(flipud(hot))
h=colorbar;
ylabel(h,'5^{o}S-5^{o}N mean extinction at 550nm (log)')
ylabel('Altitude (km a.s.l)')
xlabel('Time (year)')


%Figure 4: Scattering asymmetry factor at 550nm and 20km a.s.l.
ialt=find(alt==20);
figure(4)
contourf(time_yr,lat,squeeze(asy(:,:,ialt,4))')
colormap(flipud(hot))
h=colorbar;
ylabel(h,'Scattering asymmetry factor at 550nm and 20km a.s.l.')
ylabel('Latitude (^{o}N)')
xlabel('Time (year)')

%Et voila!
disp('Et voila!')