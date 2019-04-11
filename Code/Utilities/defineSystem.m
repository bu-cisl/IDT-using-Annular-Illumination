function [ N_obj, w_NA, x, y, u0, v0, k0, u, v, k,dx_obj] = defineSystem( lambda, NA, mag, dpix_c, NsampR, illumination_na,n_r )
%Calculates useful system parameters that define how the data is processed
%lambda: wavelength of light
%NA: Numerical Aperture
%mag: magnification
%dpix_c: size of camera pixels on sensro
%NsampR: size of the image region to be processed in [row, col] (that is, (y,x))

uMax = NA/lambda; % Maximum spatial frequency set by NA (1/um)
dxNA = (1/uMax)/2; % System resolution based on the NA (um)
dpix_m = dpix_c/mag; % Effective image pixel size on the object plane (um)
FoV = NsampR*dpix_m; % FoV in the object space

% Sample size at Fourier plane (set by the image size (FoV))
if mod(NsampR(1),2) == 1 %Row
    du(1) = 1/dpix_m./(NsampR(1)-1);
else
    du(1) = 1./FoV(1);% Sampling size at Fourier plane is always = 1/FoV
end
if mod(NsampR(2),2) == 1 %Col
    du(2) = 1/dpix_m./(NsampR(2)-1);
else
    du(2) = 1./FoV(2);
end

um_idx = uMax./du; %max spatial freq/sample size at Fourier plane = number of samples in Fourier space (?)

% Generate cutoff window by NA
n = (1:NsampR(1))-round((NsampR(1)+1)/2); %vertical, row index (y) %Indices from center of window; allows szR to be rectangle (even though this is not desired)
m = (1:NsampR(2))-round((NsampR(2)+1)/2); %horizontal, col index (x)

[mm,nn] = meshgrid(m,n); %Create grid

% assume a circular pupil function, lpf due to finite NA
w_NA = double(sqrt((nn/um_idx(1)).^2+(mm/um_idx(2)).^2)<1); %Find locations of acceptable circle in k-space

% h = fspecial('gaussian',10,5);
% w_NA = imfilter(w_NA,h);

% support of OTF (optical transfer function) is 2x of ATF(NA)
Ps_otf = double(sqrt((nn/um_idx(1)).^2+(mm/um_idx(2)).^2)<2);

phC = ones(NsampR);
% phC(ridx<0.8*um_idx&ridx>0.7*um_idx) = 0.5;
% aberration modeled by a phase function
aberration = ones(NsampR);
% aberration = exp(pi/2*1i*(exp(-(mm-20).^2/50^2-(nn+40).^2/150^2))-...
%     pi/8*1i*(exp(-(mm+40).^2/100^2-(nn-80).^2/80^2))+...
%     pi/3*1i*(exp(-(mm).^2/60^2-(nn-10).^2/30^2)));
%aberration = ones(N_m);

pupil = w_NA.*phC.*aberration;

% maxium spatial frequency achievable based on the maximum illumination
% angle from the LED array and NA of the objective
um_p = max(illumination_na(:))/lambda+uMax;
% resolution achieved after freq post-processing
dx0_p = 1/um_p/2;

NA_s = um_p*lambda;
disp(['synthetic NA is ',num2str(NA_s)]);

% assume the max spatial freq of the original object
% um_obj>um_p
% assume the # of pixels of the original object
% upsample by 2 for intensity
N_obj = round(um_p./du)*2*2; %Solving for the number of samples in freq space
%Follows rate = 2*BW where BW=um_p (max spatial freq) and rate = N_obj/du
%(num samples/sample size)

% need to enforce N_obj/Np = integer to ensure no FT artifacts
N_obj = ceil(N_obj/NsampR)*NsampR;
%%%%% Trying something out
%N_obj=NsampR; %Not normal line!
% max spatial freq of the original object
um_obj = du.*N_obj/2; 

% sampling size of the object (=pixel size of the test image)
dx_obj = 1./um_obj./2; 
% spatial coordiates for object space
[x,y] = meshgrid((-N_obj(2)/2:N_obj(2)/2-1).*dx_obj(2),(-N_obj(1)/2:N_obj(1)/2-1).*dx_obj(1));
% [xn,yn] = meshgrid(-1/2:1/N_obj(2):1/2-1/N_obj(2),-1/2:1/N_obj(1):1/2-1/N_obj(1));
%Spatial coordinates for image space
[x0,y0]=meshgrid((-NsampR(2)/2:NsampR(2)/2 - 1).*dpix_m,(-NsampR(1)/2:NsampR(1)/2 - 1).*dpix_m);
dx_img=dpix_m; %sampling size in image (demagnified pixel)

%spatial frequency coordinates
% for object space
[u,v] = meshgrid(-um_obj(2):du(2):um_obj(2)-du(2),-um_obj(1):du(1):um_obj(1)-du(1));
% for image space
[u0,v0] = meshgrid(-du(2)*NsampR(2)/2:du(2):du(2)*NsampR(2)/2-du(2),...
    -du(1)*NsampR(1)/2:du(1):du(1)*NsampR(1)/2-du(1));

%Define k-space for these coordinates
k=pi*lambda*(u.^2+v.^2)./n_r;
k0=pi*lambda*(u0.^2+v0.^2)./n_r;
xObj=(pi/lambda)*(x.^2 + y.^2); %Should all these things include n_r?
xImg=(pi/lambda)*(x0.^2 + y0.^2);

%compute depth sectioning capability
DOF0 = lambda/NA^2/2;
DOFs = lambda/NA_s^2/2;

%resolution
res0 = lambda/NA/2;
res_s = lambda/NA_s/2;

end

