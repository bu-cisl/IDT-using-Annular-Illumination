function [XY] = pol2Cart( RTh, XYmid )
%Convert polar to cartesian
%r_theta=[r, theta]; r in pixels, theta in degrees
%Output x,y coordinates in pixels

xC=RTh(:,1).*cosd(RTh(:,2))+XYmid(1);
yC=RTh(:,1).*sind(RTh(:,2))+XYmid(2);
XY=[xC yC];

end

