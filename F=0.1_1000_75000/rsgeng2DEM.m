function [f,x,y] = rsgeng2DEM(N,rL,clx,cly)
%
% [f,x,y] = rsgeng2DEM(N,rL,h,clx,cly)
%
% generates a square 2-dimensional random rough surface f(x,y) with NxN 
% surface points. The surface has a Gaussian height distribution function 
% and Gaussian autocovariance functions (in both x and y), where rL is the 
% length of the surface side, h is the RMS height and clx and cly are the 
% correlation lengths in x and y. Omitting cly makes the surface isotropic.
%
% Input:    N   - number of surface points (along square side)
%           rL  - length of surface (along square side)
%           clx, (cly)  - correlation lengths (in x and y)
%
% Output:   f  - surface heights
%           x  - surface points
%           y  - surface points
%
% Last updated: 2010-07-26 (David Bergstr�m).  
%

%format long;

x = linspace(-rL/2,rL/2,N); y = linspace(-rL/2,rL/2,N);
[X,Y] = meshgrid(x,y); 

randnsettings=rng;
save('randnsettings.mat','randnsettings')

data=dlmread('ks_all_nov.txt');
data=sort(data);
[Fi,xi] = ecdf(data);
xj = xi(2:end);
[n,m]=size(xj);
Fj = (Fi(1:end-1)+Fi(2:end))/2;
Fj = [0; Fj; 1];
xj = [0;xj;xj(n-2)+(1-Fj(n-1))*((xj(n-2)-xj(n-3))/(Fj(n)-Fj(n-1)))];
Finv = @(u) interp1(Fj,xj,u,'linear','extrap');
u = rand(N,N);

Z = Finv(u); % uncorrelated Gaussian random rough surface distribution
                   % with mean 0 and standard deviation h


    
% non-isotropic surface
if nargin == 4
    
    % Gaussian filter
    F = exp(-(X.^2/(clx^2/2)+Y.^2/(cly^2/2)));

    % correlated surface generation including convolution (faltning) and inverse
    % Fourier transform and normalizing prefactors
    f = 2/sqrt(pi)*rL/N/sqrt(clx)/sqrt(cly)*ifft2(fft2(Z).*fft2(F));
    
end
