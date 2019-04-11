function [ FIdiv, FIdivG,FI, w_2NA ] = calFI( I,xI,yI,XYmid,radP,sigmaG )
%Creates the Fourier transform images that will be processed throughout
%sigmaG generally = 2
F = @(x) fftshift(fft2(ifftshift(x))); %Fourier Transform

%Get Fourier transform images
FI=F(I);
avgFI=mean(abs(FI),3);%Average |F(I)|

w_2NA=sqrt((xI-XYmid(1)).^2 + (yI-XYmid(2)).^2)<2.*radP; %Outside 2NA support
mT=mean(avgFI(~w_2NA)); %Mean value outside 2NA support
avgFI2=max(avgFI,3.*mT); %Put a floor that will get rid of noise outside the 2NA circle
%(That is, where avgFI<3mT, make it be 3mT)

FIdiv=FI./repmat(avgFI2,[1 1 size(FI,3)]); %Divide out the average value
FIdivG=imgaussfilt(abs(FIdiv),sigmaG); %Gaussian filter it, sigma = 2

end

