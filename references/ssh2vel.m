function [ug,vg] = ssh2vel(psi,lon,lat,npts)
% SSH2VEL uses centered differences to get geostrophic velocities
% from the sea surface height fields
% 
% psi is the ssh field; psi(lat,lon)
% lon,lat are the grid, as vectors
% npts  number of points in the stencil, default is 3
% History:
% Coded by Rob Scott Oct. 2007. Confirmed that it produces geostrophic
% velocities from Aviso dynamic topography as the u,v provided by Aviso.
% Jan 2011, changed from ddc_of_2d to ddc_of_map_npts for greater accuracy
% Hardwired to use non-uniform grid but this can be changed with the last
% argument to ddc_of_map_npts.
%
% See also
% DDC_OF_MAP_NPTS

% Default is 3-pt stencil
if ~exist('npts','var')
    npts = 3;
end

if isempty(npts)
    npts = 3;
end

 geophysical_constants

% -------
% Step 1:	turn lon and lat into column vectors
% -------
 if size(lon,1) == 1
    lon = lon';
 end
 if size(lat,1) == 1
    lat = lat';
 end
 
 if size(lon,2) > 1 | size(lat,2) > 1
    error('input should be a vector for longitude latitude points')
 end 

% -------
% Step 2:	Coriolis
% -------

 [X,Y] = meshgrid(lon,lat);
 f = 2*omega*sind(Y);

% -------
% Step 3: 	Geostrophic velocities
% -------

 switch npts
     case 3
        ug = -ddc_of_map_npts(psi,lon,lat,1,npts,1)* g ./ f;
        vg =  ddc_of_map_npts(psi,lon,lat,2,npts,1)* g ./ f;
     case 5
        ug3 = -ddc_of_map_npts(psi,lon,lat,1,3,1)* g ./ f;
        vg3 =  ddc_of_map_npts(psi,lon,lat,2,3,1)* g ./ f;
        ug = -ddc_of_map_npts(psi,lon,lat,1,npts,1)* g ./ f;
        vg =  ddc_of_map_npts(psi,lon,lat,2,npts,1)* g ./ f;
        b = find(isnan(ug) & ~isnan(ug3));
        ug(b) = ug3(b);
        b = find(isnan(vg) & ~isnan(vg3));
        vg(b) = vg3(b);
     case 7
        ug3 = -ddc_of_map_npts(psi,lon,lat,1,3,1)* g ./ f;
        vg3 =  ddc_of_map_npts(psi,lon,lat,2,3,1)* g ./ f;
        ug5 = -ddc_of_map_npts(psi,lon,lat,1,5,1)* g ./ f;
        vg5 =  ddc_of_map_npts(psi,lon,lat,2,5,1)* g ./ f;
        ug = -ddc_of_map_npts(psi,lon,lat,1,npts,1)* g ./ f;
        vg =  ddc_of_map_npts(psi,lon,lat,2,npts,1)* g ./ f;
        b = find(isnan(ug) & ~isnan(ug5));
        ug(b) = ug5(b);
        b = find(isnan(vg) & ~isnan(vg5));
        vg(b) = vg5(b);
        b = find(isnan(ug) & ~isnan(ug3));
        ug(b) = ug3(b);
        b = find(isnan(vg) & ~isnan(vg3));
        vg(b) = vg3(b);
     otherwise
         error('npts must be 3, 5, or 7 or left blank so it defaults to 3')
 end
 return
 end
