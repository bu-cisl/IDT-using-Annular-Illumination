function [RTh] = cart2Pol( XY, XYmid )
%Convert polar to cartesian
%Output r_theta=[r, theta]; r in pixels, theta in degrees
%Input x_y x,y coordinates in pixels


%Theta, distance of circle center from center frequency
centD=sqrt(sum((XY-repmat(XYmid,[size(XY,1) 1])).^2,2)); %(pixels)
theta=atan2d(XY(:,2)-XYmid(2),XY(:,1)-XYmid(1)); %(degrees)

RTh=[centD, theta];

end

