% 	Written by Thomas J. Aubry, June 2019.
% 	Department of Geography, University of Cambridge
%   E-mail: ta460@cam.ac.uk
% 	Please cite the corresponding paper if you use this script

function [gmsaod, saod, reff, ext, ssa, asy, lat, alt]=postproc(SO4mass,modelpara,mstar,R_reff,wl_req)


%==========================================================================
%1) Calculate global mean SAOD and area-weighted AOD at 525nm, and
%aerosol effective radius
%==========================================================================
totmass=sum(SO4mass,2);%total mass of sulfate in the stratosphere
gmsaod525_lin=modelpara(1)*totmass;%global mean SAOD at 525nm calculated
%using the linear scaling

gmsaod525=gmsaod525_lin;
gmsaod525(totmass>mstar)=modelpara(1)*mstar^(1/3)*totmass(totmass>mstar).^(2/3);
%global mean SAOD at 525nm with 2/3 scaling applied for sulfate mass larger
%than the critical mass mstar.

waod=modelpara(1)*SO4mass.*repmat(gmsaod525./gmsaod525_lin,[1 8]);
%area-weighted AOD at 525nm in each box, calculated as the mass of sulfate
%in each box multiplied by the SAOD-sulfate mass scaling factor, and
%corrected by the ratio of the actual global mean SAOD at 525nm and the
%obtained from a linear scaling. This correction insures that the sum of
%area-weighted AOD in all boxes is equal to the global mean SAOD.

gmreff=R_reff*totmass.^(1/3);
%Calculate global mean mass-weighted effective radius in um from the
%scaling described in the companion paper.



%==========================================================================
%2) Calculate altitude and latitude dependent extinction at 525nm and
%effective radius
%==========================================================================

load('shapefunctions.mat')%load shape functions
%define latitude/altitude grid of the shape functions
lat=(-87.5:5:87.5)';
alt=(5:0.5:39.5)';

ext525=NaN(length(totmass),36,70);%pre-allocate space for extinction at 525nm
massdist=NaN(length(totmass),36,70);%pre-allocate space for sulfate mass
for i=1:length(totmass)
    %At each timestep, calculate the extinction at 525nm as the sum, over the
    %8 boxes, of the product of the area-weighted AOD in the box by the shape
    %function of the same box (cf. companion paper for more details on these
    %shape functions and how they were derived). The shape functions return
    %extinction in /km
    ext525(i,:,:)=squeeze(nansum(permute(repmat(waod(i,:)',[1 36 70]),[2 3 1]).*shapefunctions,3));
    
    %Assume that spatial distribution of sulfate mass is the same as that
    %of extinction
    massdist(i,:,:)=squeeze(nansum(permute(repmat(SO4mass(i,:)',[1 36 70]),[2 3 1]).*shapefunctions,3));
end


%Calculate weight (cos(latitude)*mass) to calculate global mean
%mass-weighted average
latweight=permute(repmat(cosd(lat),[1 70 length(totmass)]),[3 1 2]);
latweight=latweight.*massdist./repmat(sum(sum(latweight.*massdist,2),3),1,size(massdist,2),size(massdist,3));

%Assume that local effective radius follows the same spatial distribution
%as sulfate mass, raised to power 1/3.
reff=massdist.^(1/3);


%Re-scale the effective radius so that the global mean average follows the
%scaling introduced in the paper, with a minimum value of 0.1um for local
%effective radius
gmreff_scale=squeeze(nansum(nansum(reff.*latweight,2),3));
reff=0.101+reff.*repmat((gmreff-0.101)./gmreff_scale,[1 36 70]);



%==========================================================================
%3) Calculate time, altitude, latitude and wavelength dependent extinction,
%stratospheric aerosol optical depth, single scattering albedo and
%scattering asymmetry factor
%==========================================================================

%Read Mie look-up tables:
ncid = netcdf.open('eva_Mie_lookuptables.nc');

%a) Read effective radius and wavelength grid, and reformat them for
%use as input for the 2D interpolation function later
reffgrid_mie = netcdf.getVar(ncid,0);
wlgrid_mie = netcdf.getVar(ncid,1);
[Xwl,Yreff] = meshgrid(wlgrid_mie,reffgrid_mie);

%b) Read calculated parameters, which are 2D array, with one dimension for
%effective radius and the other one for wavelength.
extrat_mie = netcdf.getVar(ncid,2);%ratio of extinction to extinction at 550nm (EXT)
ssa_mie = netcdf.getVar(ncid,3);%single scattering albedo (SSA)
asy_mie = netcdf.getVar(ncid,4);%scattering asymmetry factor (ASY)

%preallocate memory for calculating EXT, SSA and ASY
ext=NaN(size(ext525,1),size(ext525,2),size(ext525,3),length(wl_req));
ssa=NaN(size(ext));
asy=NaN(size(ext));

%c) Loop through latitude, altitude and wavelength to calculate them. All
%calculations are done by linearly interpolating the Mie lookup tables at
%the requested wavelength and the effective radius outputted by the model.
for ilat=1:length(lat)
    for ialt=1:length(alt)
        for iwl=1:length(wl_req)
            %ignore points in time where the extinction at 525nm or
            %effective radius are NaNs
            mask=~isnan(squeeze(ext525(:,ilat,ialt))) & ~isnan(squeeze(reff(:,ilat,ialt)));
            %Calculate the raio of extinction at desired wavelength to
            %extinction at 525nm
            ratio525=interp2(Xwl,Yreff,extrat_mie,wl_req(iwl),min(reff(mask,ilat,ialt),1.29))./interp2(Xwl,Yreff,extrat_mie,0.525,min(reff(mask,ilat,ialt),1.29));
            %multiple above ratio by extinction at 525nm
            ext(mask,ilat,ialt,iwl)=ext525(mask,ilat,ialt).*ratio525;
            %calculate SSA and ASY
            ssa(mask,ilat,ialt,iwl)=interp2(Xwl,Yreff,ssa_mie,wl_req(iwl),min(reff(mask,ilat,ialt),1.29));
            asy(mask,ilat,ialt,iwl)=interp2(Xwl,Yreff,asy_mie,wl_req(iwl),min(reff(mask,ilat,ialt),1.29));
            
            
        end
    end
end

%Calculate stratospheric aerosol optical depth. These are simply the sum of
%extinction along vertical dimension, multiplied by 0.5 because the
%vertical grid is regurlarly spaced by 0.5km. All tropospheric values are
%NaNs in the shape functions and thus in ext.
saod=squeeze(nansum(ext,3)*0.5);

%Calculate weights (cosinus(latitude)) for calculating global mean 
latweight=permute(repmat(cosd(lat),[1 length(wl_req) length(totmass)]),[3 1 2]);
latweight=latweight/squeeze(sum(latweight(1,:,1),2));
%Calculate global mean SAOD
gmsaod=squeeze(nansum(saod.*latweight,2));

end