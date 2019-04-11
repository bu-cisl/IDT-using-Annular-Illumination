function [newx, newy, newvalue] = ada_nonmax_suppression(xp, yp, value, n)
% Adaptive non-maximun suppression 
% For each Harris Corner point, the minimum suppression radius is the
% minimum distance from that point to a different point with a higher 
% corner strength. 
% Input:
% xp,yp - coordinates of harris corner points
% value - strength of suppression
% n - number of interesting points
% Output:
% newx, newy - new x and y coordinates after adaptive non-maximun suppression
% value - strength of suppression after adaptive non-maximun suppression

% ALLOCATE MEMORY
% newx = zeros(n,1);
% newy = zeros(n,1);
% newvalue = zeros(n,1);

if(length(xp) < n)
newx = xp;
newy = yp;
newvalue = value;
return;
end

radius = zeros(n,1);
c = .9;
maxvalue = max(value)*c;
for i=1:length(xp)
if(value(i)>maxvalue)
radius(i) = 99999999;
continue;
else
dist = (xp-xp(i)).^2 + (yp-yp(i)).^2;
dist((value*c) < value(i)) = [];
radius(i) = sqrt(min(dist));
end
end

[~, index] = sort(radius,'descend');
index = index(1:n);

newx = xp(index);
newy = yp(index);
newvalue = value(index);
