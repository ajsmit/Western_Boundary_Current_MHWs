function coeffs = get_npts_stencil_coeffs(xaxis,Npts)
% GET_NPTS_STENCIL_COEFFS computes the coefficients needed for Npts stencil
%
% xaxis     vector of either latitude or longitude points (not necessarily
%           x-direction)
% Npts    Number of points in stencil, must be {3,5,7}.
%
% Coefficients, cn, needed for the linear combination of Npts psi points to get
% the first derivative. Does not take into account cos(lat) when xaxis is
% longitude. This will be taken into account in DDC_OF_MAP_NPTS.
%
% History:
% Coded by Rob Scott Jan 2011. 
% References:
% Arbic et al. 2011 JGR 
%
% see also:
% DDC_OF_MAP_NPTS

geophysical_constants % need mperdeg

npts_left = (Npts-1)/2;   % = 1 for Npts = 3, etc.
npts_right = npts_left;

if size(xaxis,1) == 1
    xaxis = xaxis';
end

if size(xaxis,2) > 1
    error('input should be a vector for latitude points')
end

% Initialize storage varibles
Ny     = length(xaxis);
coeffs = NaN(Ny,Npts);

switch Npts
    case 3   %      
        for jj = 1+npts_left:Ny-npts_left
            A    = [ones(1,Npts) ; ...
                xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj); ...
                (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^2];
            B    = zeros(Npts,1);
            B(2) = 1;
            tmp  = A\B;
            coeffs(jj,:) = tmp'/mperdeg;
        end
    case 5
        for jj = 1+npts_left:Ny-npts_left
            A    = [ones(1,Npts) ; ...
                    xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj); ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^2; ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^3; ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^4];
            B    = zeros(Npts,1);
            B(2) = 1;
            tmp  = A\B;
            coeffs(jj,:) = tmp'/mperdeg;
        end        
    case 7
        for jj = 1+npts_left:Ny-npts_left
            A    = [ones(1,Npts) ; ...
                    xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj); ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^2; ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^3; ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^4; ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^5; ...
                    (xaxis(jj-npts_left:jj+npts_left)'-xaxis(jj)).^6];
            B    = zeros(Npts,1);
            B(2) = 1;
            tmp  = A\B;
            coeffs(jj,:) = tmp'/mperdeg;
        end
    otherwise
        error('Npts must be 3,5, or 7')
end