function [ freqXY, con, radP, xI, yI, uI, vI, XYmid ] = calCoord(freqUV,imSz,dpix_c,mag,NA,lambda )
%Convert from k-space coordinates to pixel coordinates
con=imSz(1).*dpix_c./mag; %Conversion factor (pixels/(1/um))

%k-space (u,v) coordinates (1/um)
uCent=freqUV(:,1); 
vCent=freqUV(:,2);

%Real space (x,y) coordinates (pixels)
xMid=floor(imSz(2)./2)+1; 
yMid=floor(imSz(1)./2)+1; 

[xI,yI]=meshgrid(1:imSz(2),1:imSz(1)); %Grid

xCent=xMid+uCent.*con; %k-space coordinates in terms of pixel values
yCent=yMid+vCent.*con;

%Combine into variables for ease of transport
freqXY=[xCent, yCent];
XYmid=[xMid, yMid];

%Theta, distance of circle center from center frequency
%freqDTh=cart2Pol(freqXY, XYmid);

%Predicted radius (pixels)
radP=(NA/lambda).*con; 

%k-space (u,v) coordinates (1/um)
uI=(xI-xMid)./con;
vI=(yI-yMid)./con;


end

