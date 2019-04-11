function [hh, inliers] = ransacfithomography(ref_P, dst_P, npoints, threshold);
% 4-point RANSAC fitting
% Input:
% matcher_A - match points from image A, a matrix of 3xN, the third row is 1
% matcher_B - match points from image B, a matrix of 3xN, the third row is 1
% thd - distance threshold
% npoints - number of samples
% 
% 1. Randomly select minimal subset of points
% 2. Hypothesize a model
% 3. Computer error function
% 4. Select points consistent with model
% 5. Repeat hypothesize-and-verify loop
% 
% Yihua Zhao 02-01-2014
% zhyh8341@gmail.com

ninlier = 0;
fpoints = 8; %number of fitting points
for i=1:2000
rd = randi([1 npoints],1,fpoints);
pR = ref_P(:,rd);
pD = dst_P(:,rd);
h = getHomographyMatrix(pR,pD,fpoints);
rref_P = h*ref_P;
rref_P(1,:) = rref_P(1,:)./rref_P(3,:);
rref_P(2,:) = rref_P(2,:)./rref_P(3,:);
error = (rref_P(1,:) - dst_P(1,:)).^2 + (rref_P(2,:) - dst_P(2,:)).^2;
n = nnz(error<threshold);
if(n >= npoints*.95)
hh=h;
inliers = find(error<threshold);
pause();
break;
elseif(n>ninlier)
ninlier = n;
hh=h;
inliers = find(error<threshold);
end 
end
