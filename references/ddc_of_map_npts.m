function deriv = ddc_of_map_npts(psi,lon,lat,div_der,Npts,flag_grid)
% DDC_OF_MAP_NPTS partial derivative a 2d mapped field using Npts stencil
% 
% Calculates the first derivative of psi with respect to either x or y
% using 3,5, or 7 point stencils. It assumes that psi is a mapped field on
% the globe, so it is periodic in the x-direction. It exploits this
% peridicity so that there are not missing data points at the dateline or
% Greenwich meridians. 
%
% Input:
% psi is the field; psi(lat,lon)
% lon,lat are the grid, as vectors
% div_der == 1, \partial psi/  partial lat
% div_der == 2, \partial psi/  partial lon
% Npts    Number of points in stencil, must be {3,5,7}.
% flag_grid  = 0 uniform; = 1 non-uniform
%
% Output:
% deriv = \frac{\partial psi}{\partial x_{div_der}} where
% x_1 = y 
% x_2 = x
% Units of deriv = units of psi/[m]
%
% History:
% Coded by Rob Scott Jan 2011. 
% References:
% Arbic et al. 2011 JGR 
%
% see also
% GET_NPTS_STENCIL_COEFFS
% TEST_NPT_DER

geophysical_constants

% -------
% Step 1: turn lon and lat into column vectors
% -------
% 
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
% Step 2: Stencil parameters:
% -------
% 

npts_left = (Npts-1)/2;   % = 1 for Npts=3, etc.
npts_right = npts_left;

% -------
% Step 3: (for d/dx only) since it's periodic in x
% -------
% make cyclic, padding with width sufficient for the stencile used.
 if div_der == 2
    nx_in = length(lon);
    lon = [ lon(end-npts_left+1:end)-360; lon; 360+lon(1:npts_right)];
    psi = [psi(:,end-npts_left+1:end) psi psi(:,1:npts_right)];
 end
 
 % -------
 % Step 4: Get coefficients of the Npt-stencil
 % -------
 % 
 
 if Npts > 0
     if div_der == 1
         coeffs = get_npts_stencil_coeffs(lat,Npts);
     elseif div_der == 2
         coeffs = get_npts_stencil_coeffs(lon,Npts);
     end
 end
 
 % -------
 % Step 5: Get derivative depending on direction and grid and stencil
 % -------
 
 [X,Y] = meshgrid(lon,lat);
 
 [Ny,Nx] = size(psi);
 
 % Initialize storage varibles
 
 deriv  = NaN(Ny,Nx);
 
 if div_der == 1
     if flag_grid == 0 % for historical reasons, keep this method
         switch Npts
             case 3
                 deriv =      (psi(3:Ny,:) - psi(1:Ny-2,:))  ...
                     ./ (mperdeg*(Y(3:Ny,:) -   Y(1:Ny-2,:)));
             case 5
                 deriv =(-psi(5:Ny,:)+8*psi(4:Ny-1,:)-8*psi(2:Ny-3,:)+psi(1:Ny-4,:)) ...
                     ./ (mperdeg*(Y(4:Ny-1,:) - Y(2:Ny-3,:))*6);
             otherwise
                 error('only 3 pt or 5 pt stenciles allowed for uniform grid')
         end
% pad with NaN so size(deriv) = size(psi)
         deriv = [NaN(npts_left,Nx); deriv; NaN(npts_right,Nx)];
     elseif flag_grid == 1 % preferred method, account for non-uniform grid
         for jj = 1+npts_left:Ny-npts_left
             %coeff_band: 1x Npts
             coeff_band = coeffs(jj,:);
             % psi_band: Npts x Nx
             psi_band   = psi(jj-npts_left:jj+npts_left,:);
             % deriv_band: 1 x Nx
             deriv_band =  coeff_band * psi_band ;
             deriv(jj,:) = deriv_band;
         end
         
     end
 elseif div_der ==2
     if flag_grid == 0 % for historical reasons, keep this method
         switch Npts
             case 3
                 % psi_x
                 deriv =      (psi(:,3:Nx) - psi(:,1:Nx-2)) ...
                     ./ (mperdeg*(X(:,3:Nx) -   X(:,1:Nx-2)) .* cosd(Y(:,2:Nx-1)) );
             case 5
                 deriv =(-psi(:,5:Nx)+8*psi(:,4:Nx-1)-8*psi(:,2:Nx-3)+psi(:,1:Nx-4)) ...
                     ./ (mperdeg*(X(:,4:Nx-1)-X(:,2:Nx-3)).*cosd(Y(:,3:Nx-2))*6);
             otherwise
                 error('only 3 pt or 5 pt stenciles allowed for uniform grid')
         end
     elseif flag_grid == 1 % lon is uniform, but for debugging also use this general
         % nonuniform code -- it should give same results as flag_grid == 0
         for ii = 1+npts_left:Nx-npts_left
             %coeff_band: Npts x 1;
             coeff_band = coeffs(ii,:)';
             % psi_band: Ny x Npts
             psi_band   = psi(:,ii-npts_left:ii+npts_left);
             % deriv_band: Ny x 1
             deriv_band =  psi_band * coeff_band;
             deriv_band =  deriv_band ./ cosd(lat);  % take cos lat factor into account
             deriv(:,ii) = deriv_band;
         end
         
         % trim redundancy in x-dir
         deriv = deriv(:,1+npts_left:end-npts_right);
         
     end
 else
     disp('4th arg. must be 1 or 2 for y or x direction')
 end

