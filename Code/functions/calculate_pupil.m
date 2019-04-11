function conds=calculate_pupil(NA,pixelsize,k,M,N)


n=1:N;
m=1:M;
L0X=pixelsize*M;
L0Y=pixelsize*N;

v=(-M/L0X/2+1/L0X*(m-1));
u=(-N/L0Y/2+1/L0Y*(n-1));

[uu,vv] = meshgrid(u,v); 
SH=uu.*uu+vv.*vv;


kk=k/(2*pi)*NA;
conds=zeros(M,N);
conds(SH<=kk^2)=1;

end