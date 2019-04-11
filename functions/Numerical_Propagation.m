%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Uz  = Numerical_Propagation(U0,z,pixelsize,lambda,NA,method)
% function  Uz  = Numerical_Propagation(U0,z,pixelsize,lambda,method)
% Purpose: numerical propagation the complex field to another plane at a
% given distance using 'Angular Specturm' or 'Fresnel' method
%'Inputs': 'U0',Orignal complex field;
%          'z',Propagation distance 
%          'pixelsize'£¬pixelsize('mm'); 
%          'lambda',Wavelength('mm'); 
%          'method', type of transfer function used ('Angular Spectrum' or
%          'Fresnel')
%'Outputs': 'Uz',Complex field after propagation; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 2.0 - 
% Coded by Chao Zuo - 2012-11-9   
% Lastest edited by Chao Zuo - 2014-7-15 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(U0);

x=1:M;
y=1:N;

L0X=pixelsize*M;
L0Y=pixelsize*N;

k=2*pi/lambda;

u=lambda*(-M/L0X/2+1/L0X*(x-1));
v=lambda*(-N/L0Y/2+1/L0Y*(y-1));

[uu,vv] = meshgrid(u,v); 


SHy=(-(M)/2:(M)/2-1);
SHy=SHy./(M*pixelsize);
SHy=ones(N,1)*SHy;
SHy=SHy.*SHy;

SHx=(-(N)/2:(N)/2-1);
SHx=SHx./(N*pixelsize);
SHx=SHx'*ones(1,M);
SHx=SHx.*SHx;
SH=SHx+SHy;

kk=k/(2*pi)*NA;
HHH=zeros(N,M);
HHH(SH<=kk^2)=1;

FU0=fftshift(fft2(U0));

if(strcmp(method,'Angular Spectrum'))
H=exp(1i*k*z*sqrt(1-uu.^2-vv.^2)); % Angular Specturm method 
elseif(strcmp(method,'Fresnel'))
H=exp(1i*k*z*(1-(uu.^2+vv.^2)/2)); % Fresnel method
else
errordlg('Type of transfer function must be <Angular Spectrum> or <Fresnel>','Error');
end

Uz=ifft2(fftshift(FU0.*H.*HHH));


%%
