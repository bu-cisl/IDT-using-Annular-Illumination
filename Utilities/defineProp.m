function [ H0 ] = defineProp( lambda, z0, k0, mode)
%Propagates by specified mode
%Used if sample is not precisely in focus

F = @(x) ifftshift(fft2(fftshift(x))); %Fourier Transform

% Fresnel
% H2 = exp(1i*2*pi/lambda*dz)*exp(-1i*pi*lambda*dz*(u.^2+v.^2));
% angular spectrum
% H = exp(1i*2*pi*sqrt(1/lambda^2-u.^2-v.^2)*dz);
% mask_evanescent = double(sqrt(u0.^2+v0.^2)<1/lambda);
% H = H.*mask_evanescent;

if mode == 0 %Fresnel
    %Make sure we don't alias
    %H0 = exp(1i*z0*2*pi*1/lambda).*exp(-1i*k0*z0); 
    H0 =exp(-1i*k0*z0);  %Added refractive index?


else %Otherwise, angular spectrum
    eva0 = double(k0<pi/lambda);
    H0 = exp(1i*2*pi*sqrt((1/lambda^2-k0/pi/lambda).*eva0)*z0);
    
end


end

