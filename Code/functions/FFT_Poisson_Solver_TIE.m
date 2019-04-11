%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Solution to TIE using FFT-based Possion solver
%'Inputs': 'dIdz',Intensity derivative along optical axis;
%          'I0', Infoucs intensity image 
%          'pixelsize'£¬pixelsize('mm'); 
%          'k',Wavenumber
%          'r',regularzation parameter (to remove low-frequency noise)
%'Outputs': 'phi',object phase map; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = FFT_Poisson_Solver_TIE(dIdz,I0,pixelsize,k,r)

[M,N]=size(dIdz);
n=1:N;
m=1:M;
L0X=pixelsize*M;
L0Y=pixelsize*N;

v=(-M/L0X/2+1/L0X*(m-1));
u=(-N/L0Y/2+1/L0Y*(n-1));

[uu,vv] = meshgrid(u,v); 

kdIdz=-k*dIdz;

Fleft=fft2(kdIdz);

Fphi=fftshift(Fleft).*(-4*pi*pi*(uu.*uu+vv.*vv))./(r+(-4*pi*pi*(uu.*uu+vv.*vv)).^2);
bigphi=real(ifft2(fftshift(Fphi)));
Fbigphi=fft2(bigphi);

Fphi=fftshift(Fbigphi).*(2*1i*pi*(uu));
dxbigphi=real(ifft2(fftshift(Fphi)));


Fphi=fftshift(Fbigphi).*(2*1i*pi*(vv));
dybigphi=real(ifft2(fftshift(Fphi)));

dxbigphi=dxbigphi./I0;
dybigphi=dybigphi./I0;

Fbigphi=fft2(dxbigphi);

Fphi=fftshift(Fbigphi).*(2*1i*pi*(uu));
dxdxbigphi=real(ifft2(fftshift(Fphi)));


Fbigphi=fft2(dybigphi);

Fphi=fftshift(Fbigphi).*(2*1i*pi*(vv));
dydybigphi=real(ifft2(fftshift(Fphi)));

ddphi=dxdxbigphi+dydybigphi;

Fleft=fft2(ddphi);

Fphi=fftshift(Fleft).*(-4*pi*pi*(uu.*uu+vv.*vv))./(r+(-4*pi*pi*(uu.*uu+vv.*vv)).^2);

phi=real(ifft2(fftshift(Fphi)));



