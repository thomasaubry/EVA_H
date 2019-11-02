% 	Written by Thomas J. Aubry, June 2019.
% 	Department of Geography, University of Cambridge
%   E-mail: ta460@cam.ac.uk
% 	Please cite the corresponding paper if you use this script


%Note: the structure and indices of model boxes to which we refer in this
%script are defined in the companion paper. Boxes are indexed from 1 to 8,
%from South to North and from top to bottom. 

%Given the boundaries of the boxes of our model (defined by h1lim, h2lim
%and latlim) and a path to an excel file containing a list of volcanic SO2
%emission parameters, this function returns the mass of sulfur injected in
%each box of the model and the time of injection.


function [injec timelist]=so2injection_8boxes(h1lim,h2lim,latlim,filepath)
%injec will be a N x 8 matrix containing the mass of sulfur to be injected
%in each box in Tg S (N=number of eruptive events).
%timelist will contain the time of injection, in months after Jan 1st 1979
%filepath must be a relative path to a .xlsx file where column and units
%are the same as the two example files provided (singleinjection.xlsx and
%volSO2_Carn2016.xlsx)

%==========================================================================
%Load data for tropopause height from NCEP/NCAR reanalysis
%==========================================================================

load('ncep_tropo.mat'); %zonal mean tropopause height for 1979-2016 
%Define and time corresponding to NCEP tropopause
lat=-90:2.5:90;
tropotime=0.5:1:455.5; %time in month since 1979
%Define average tropopause height in latitudinal band of the model
[a ntl]=min(abs(lat-latlim));
[a stl]=min(abs(lat+latlim));
climtropoheight=NaN(size(tropoheight,2),3);
climtropoheight(:,1)=sum(tropoheight(1:stl-1,:)'.*(repmat(cosd(lat(1:stl-1)),[size(tropoheight,2) 1])/sum(cosd(lat(1:stl-1)))),2);
climtropoheight(:,2)=sum(tropoheight(stl:ntl,:)'.*(repmat(cosd(lat(stl:ntl)),[size(tropoheight,2) 1])/sum(cosd(lat(stl:ntl)))),2);
climtropoheight(:,3)=sum(tropoheight(ntl+1:end,:)'.*(repmat(cosd(lat(ntl+1:end)),[size(tropoheight,2) 1])/sum(cosd(lat(ntl+1:end)))),2);

%==========================================================================
%Load SO2 injection parameters
%==========================================================================
[filedata filetext]=xlsread(filepath);
so2mass=filedata(:,6);%SO2 mass in kt of SO2
erulat=filedata(:,4);%eruption latitude in degree N
eruhstar=(filedata(:,8)./filedata(:,7));%ratio of injection height/tropopause height
eruheight=(filedata(:,8));%injection height in km a.s.l.
erutropo=filedata(:,7);%tropopause height in km a.s.l.
erudate=filedata(:,1:3);%date
erudate=(datenum(erudate)-datenum([1979 1 1 0 0 0]))*12/365.25; %date in number of months since 1979
timelist=erudate;
so2mass=(10^(-3))*(0.50052)*so2mass; %convert mass from kt so2 to tg of S


%==========================================================================
%Main loop
%==========================================================================

%Pre-allocate memory for list of mass injected in each box
injec=zeros(length(erudate),8);

%Loop through all eruptions
for i=1:length(erudate)

%==========================================================================
%Define parameters and distribution function for SO2 injection #i
%==========================================================================  
h0=eruheight(i);%height
l0=erulat(i);%latitude
dh=1.2;%vertical thickness of the eruption cloud, in km
dl=7;%latitudinal extent of the eruption cloud, in degree
so2latdist = @(l) exp(-((l-l0)./dl).^2);%vertical distribution function (gaussian)
so2hgtdist = @(h) exp(-((h-h0)./dh).^2);%latitudinal distribution function  (gaussian)

%Find height of tropopause in latitudinal band
[a ert]=min(abs(erudate(i)-tropotime));
SHtrop=climtropoheight(ert,1);
Ttrop=climtropoheight(ert,2);
NHtrop=climtropoheight(ert,3);

%For the band where the eruption occur, use local tropopause height at
%volcano location
if erulat(i)<-latlim
    SHtrop=erutropo(i);
elseif erulat(i)>latlim
 NHtrop=erutropo(i);
else
 Ttrop=erutropo(i);   
end


%==========================================================================
%For each box, calculate the mass injected by eruption #i
%========================================================================== 


%Box 1

%1) set integration limits for the box
lmin=-inf; lmax=-latlim; %latitudinal limits; box 1 comprise latitudes south of -latlim
hmin=h2lim; hmax=h0+100*dh; %vertical limits; box one is above h2lim

%2) Calculate injected mass in this box by multiplying the mass injected by
%the eruption by the product of the distribution functions integrated over the box
%limit, and normalized by their integrals over the entire domain.
%Essentially, this is calculating the fraction of mass injected in a box
%and multiplying by the mass injected by the eruption.
injec(i,1) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));


%Box 2
lmin=-latlim; lmax=+latlim;
hmin=h2lim; hmax=h0+100*dh;
injec(i,2) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));

%Box 3
lmin=latlim; lmax=l0+100*dl;
hmin=h2lim; hmax=h0+100*dh;
injec(i,3) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));


%Box 4
lmin=l0-100*dl; lmax=-latlim;
hmin=Ttrop; hmax=h2lim; %the bottom of boxes 4-6 is the tropical tropopause height
injec(i,4) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));

%Box 5
lmin=-latlim; lmax=+latlim;
hmin=Ttrop; hmax=h2lim;
injec(i,5) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));

%Box 6
lmin=latlim; lmax=l0+100*dl;
hmin=Ttrop; hmax=h2lim;
injec(i,6) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));

%Box 7
lmin=-l0-100*dl; lmax=-latlim;
hmin=SHtrop; hmax=Ttrop;%the bottom of box 7 is the SH tropopause height
injec(i,7) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));


%Box 8
lmin=latlim; lmax=l0+100*dl;
hmin=NHtrop; hmax=Ttrop; %the bottom of box 8 is the NH tropopause height
injec(i,8) = so2mass(i)*(integral(so2latdist,lmin,lmax)*integral(so2hgtdist,hmin,hmax))/(integral(so2hgtdist,h0-100*dh,h0+100*dh)*integral(so2latdist,l0-100*dl,l0+100*dl));

end

injec=injec';

end