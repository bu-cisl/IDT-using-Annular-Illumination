function [descriptors] = getFeatureDescriptor(input_image, xp, yp, sigma)
% Extract non-rotation invariant feature descriptors
% Input:
% input_image - input gray-scale image
% xx - x coordinates of potential feature points
% yy - y coordinates of potential feature points
% output:
% descriptors - array of descriptors

% FIRST BLUR WITH GAUSSIAN KERNEL
g = fspecial('gaussian', 5, sigma);
blurred_image = imfilter(input_image, g, 'replicate','same');

% THEN TAKE A 40x40 PIXEL WINDOW AND DOWNSAMPLE TO 8x8 PATCH
npoints = length(xp);
descriptors = zeros(npoints,64);

for i = 1:npoints
%pA = imresize( blurred_image(xp(i)-20:xp(i)+19, yp(i)-20:yp(i)+19), .2);
patch = blurred_image(xp(i)-20:xp(i)+19, yp(i)-20:yp(i)+19);
patch = imresize(patch, .2);
descriptors(i,:) = reshape((patch - mean2(patch))./std2(patch), 1, 64); 
end
