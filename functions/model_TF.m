function [PTF_Output,ATF_Output] = model_TF(source_pattern,deltaz,Pixelsize,lambda,NA)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Nx,Ny]=size(source_pattern);

AA_Pha=ones(Nx,Ny);
PP_Pha=zeros(Nx,Ny);
PP_Pha(Nx/2+1,Ny/2+1)=0.1;
PurePha_Obj=AA_Pha.*exp(1i*PP_Pha);

AA_Abs=ones(Nx,Ny);
AA_Abs(Nx/2+1,Ny/2+1)=0.9;
PP_Abs=zeros(Nx,Ny);
PureAbs_Obj=AA_Abs.*exp(1i*PP_Abs);

S=ifft2(fftshift(source_pattern));
SS=angle(S);
U0_Pha=PurePha_Obj.*exp(1i.*SS);
U0_Abs=PureAbs_Obj.*exp(1i.*SS);

Uz1  = Numerical_Propagation(U0_Pha,deltaz,Pixelsize,lambda,NA,'Angular Spectrum');
Uz2  = Numerical_Propagation(U0_Abs,deltaz,Pixelsize,lambda,NA,'Angular Spectrum');

Iz1=abs(Uz1).^2;
Iz2=abs(Uz2).^2;

PTF_Output=single(fft2(Iz1-1)./(fft2(PP_Pha)))./2;
ATF_Output=single(fft2(Iz2+1)./(fft2(AA_Abs)))./2; 
end

