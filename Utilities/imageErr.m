function [ rmse ] = imageErr( img, xC, yC, rad, xI, yI)
%Calculates rmse of image filtered by a pupil at centers xC, yC, with
%radius rad on grid xI, yI

F = @(x) fftshift(fft2(ifftshift(x))); %Fourier Transform
Ft = @(x) fftshift(ifft2(ifftshift(x))); %Inverse Fourier Transform

%Check if need to expand matrices
if size(img,3)>numel(xC)
    if numel(xC) == 1
        xC=xC.*ones([size(img,3) 1]);
        yC=yC.*ones([size(img,3) 1]);
    else
        error('Size (xC,yC) ~= size(img,3), but xC, yC are vectors. Unclear how to proceed.')
    end
elseif size(img,3)<numel(xC)
    if size(img,3) == 1
        img = repmat(img,[1 1 numel(xC)]);
    else
        error('Size(img,3) ~= size (xC,yC), but img is 3D. Unclear how to proceed.')
    end
end

pupil=zeros(size(img));
for ii=1:size(img,3)
    pupil(:,:,ii)=sqrt((xC(ii)-xI).^2+(yC(ii)-yI).^2)<=rad;
end

Iest=Ft(F(img).*pupil);

rmse=squeeze(sqrt(sum(sum((abs(Iest)-img).^2))./numel(img(:,:,1))));

end
