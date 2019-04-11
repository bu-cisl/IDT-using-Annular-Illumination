%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Background Phase Removal
%'Inputs': 'phi_in',Input Phase;
%          'Th', Background threshond
%          'sampling', Fitting data sampling
%          'degree', Degrees of freedom
%'Outputs': 'phi', Phase map with background removed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi,back] = Remove_Background(phi_in,Th,sampling,degree)

[Nx,Ny] = size(phi_in);
[xx,yy] = meshgrid((-fix(Nx/2):1:fix((Nx-1)/2)),(-fix(Ny/2):1:fix((Ny-1)/2)));

x_data = xx(phi_in<Th);
y_data = yy(phi_in<Th);
phi_data = phi_in(phi_in<Th);

x_data = x_data(:);
y_data = y_data(:);
phi_data = phi_data(:);

sp=polyfitn([x_data(1:sampling:end),y_data(1:sampling:end)],phi_data(1:sampling:end),degree);
back=polyvaln(sp,[xx(:),yy(:)]);
back=reshape(back,Nx,Ny);
back=back-mean2(back);
phi = phi_in - back;
