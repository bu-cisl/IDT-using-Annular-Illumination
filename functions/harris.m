function [xp, yp, value] = harris(input_image, sigma,thd, r)
% Detect harris corner 
% Input:
% sigma - standard deviation of smoothing Gaussian
% r - radius of region considered in non-maximal suppression
% Output:
% xp - x coordinates of harris corner points
% yp - y coordinates of harris corner points
% value - values of R at harris corner points

% CONVERT RGB IMAGE TO GRAY-SCALE, AND BLUR WITH G1 KERNEL
g1 = fspecial('gaussian', 7, 1);
gray_image = imfilter(input_image, g1);

% FILTER INPUT IMAGE WITH SOBEL KERNEL TO GET GRADIENT ON X AND Y
% ORIENTATION RESPECTIVELY
h = fspecial('sobel');
Ix = imfilter(gray_image,h,'replicate','same');
Iy = imfilter(gray_image,h','replicate','same');

% GENERATE GAUSSIAN FILTER OF SIZE 6*SIGMA (± 3SIGMA) AND OF MINIMUM SIZE 1x1
g = fspecial('gaussian',fix(6*sigma), sigma);

Ix2 = imfilter(Ix.^2, g, 'same').*(sigma^2); 
Iy2 = imfilter(Iy.^2, g, 'same').*(sigma^2);
Ixy = imfilter(Ix.*Iy, g, 'same').*(sigma^2);

% HARRIS CORNER MEASURE
R = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps); 
% ANOTHER MEASUREMENT, USUALLY k IS BETWEEN 0.04 ~ 0.06
% response = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2;

% GET RID OF CORNERS WHICH IS CLOSE TO BORDER
R([1:20, end-20:end], :) = 0;
R(:,[1:20,end-20:end]) = 0;

% SUPRESS NON-MAX 
d = 2*r+1; 
localmax = ordfilt2(R,d^2,true(d)); 
R = R.*(and(R==localmax, R>thd));

% RETURN X AND Y COORDINATES 
[xp,yp,value] = find(R);
