function [ pixMean, pixVal ] = radialAverage( img, radius, radiusV, thetaV, distV )
%Takes the radial average from a center point defined by its theta,
%distance from the center of the image
%radius - expected radius of circle
%radiusV - vector of radii to sample in pixels
%thetaV,distV - vectors that define the circle centers as an angle (degrees) and
%distance (pixels) from the center of the image (polar coordinates)

% radV=(2:-0.5:-2)' + rad; %Radii to test
% chTh=-5:0.25:5; %Angle of circle center from origin in degrees
% chDi=-5:0.5:5; %Distance of circle center from origin
if size(radiusV,1)==1
    radiusV=radiusV';
end
if size(thetaV,1)~=1
    thetaV=thetaV';
end
if size(distV,1)~=1
    distV=distV';
end

thetaV=wrapTo180(thetaV);
[xI,yI]=meshgrid(1:size(img,2),1:size(img,1));
xMid=ceil(size(img,2)/2); 
yMid=ceil(size(img,1)/2);

numTh=length(thetaV);
numD=length(distV);
numR=length(radiusV);

%Set up angles defining circle
angleV=wrapTo180(0:180/64:360); %(in degrees), row vector
numA=length(angleV);

%Pre expand some of the matrices that don't depend on image
rMat=repmat(radiusV,[1 numA numTh numD]);
aMat=repmat(angleV,[numR 1 numTh numD]);

rcos=rMat.*cosd(aMat);
rsin=rMat.*sind(aMat);

[arcMat, numA2]=calArc(radius, distV, thetaV, angleV); %Get arc
arcMat=repmat(arcMat,[numR 1 1 numD]);

%Calculate xC, yC along axial spokes to try
xC=repmat(distV,[numTh 1]).*cosd(repmat(thetaV',[1 numD]))+xMid;
yC=repmat(distV,[numTh 1]).*sind(repmat(thetaV',[1 numD]))+yMid;

%Rad in first dim, angle in 2nd
%Put theta in 3rd dim, distance in 4th
%Calculate the values of the circles
xR=reshape(rcos(arcMat),[numR numA2 numTh numD])+repmat(permute(xC,[3 4 1 2]),[numR numA2 1 1]);
yR=reshape(rsin(arcMat),[numR numA2 numTh numD])+repmat(permute(yC,[3 4 1 2]),[numR numA2 1 1]);

%Interpolate the circle values
pixVal=interp2(xI,yI,img,xR,yR,'spline');

%Derivatives
pixMean=squeeze(mean(pixVal,2)); %Mean value around rim


end

