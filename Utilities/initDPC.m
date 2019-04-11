function [ objR ] = initDPC( I, source_list_na, illumNA, NA_obj, lambda, dpix_sample, NsampR, nObj)
%DPC initialization for FPM

F = @(x) ifftshift(fft2(fftshift(x))); %Fourier Transform
Ft = @(x) ifftshift(ifft2(fftshift(x))); %Inverse Fourier Transform
synthDPC = zeros([NsampR 4]); %Init

%Identify which quadrant each thing belongs in
BFI=illumNA<NA_obj;
synthDPC(:,:,1) = sum(I(:,:,(source_list_na(:,2)>0 & BFI)),3); %Top
synthDPC(:,:,2) = sum(I(:,:,(source_list_na(:,2)<=0 & BFI)),3); %Bottom
synthDPC(:,:,3) = sum(I(:,:,(source_list_na(:,1)>0 & BFI)),3); %Left
synthDPC(:,:,4) = sum(I(:,:,(source_list_na(:,1)<=0 & BFI)),3); %Right

[~,dpc_phase] = dpc_phaseimaging(synthDPC,NA_obj,lambda,[0,180,90,270],dpix_sample,1,1,1e-2);
ObjF = padarray(F(sqrt(mean(I(:,:,illumNA<NA_obj),3)).*exp(1i*dpc_phase)),[(nObj(1)-NsampR(1))/2,(nObj(2)-NsampR(2))/2]);
objR=Ft(ObjF);
end

