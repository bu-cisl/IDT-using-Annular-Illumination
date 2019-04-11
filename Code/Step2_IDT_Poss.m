
%%%%% IDT poss after LED Postion calibtion %%%%%

% apply calibrated
 if Calib==1
    Ini_NAx=(freqXY3(:,2)'-(Cablib_Nx/2+1))/Cablib_Nx/Pixelsize;
    Ini_NAy=(freqXY3(:,1)'-(Cablib_Ny/2+1))/Cablib_Ny/Pixelsize;
 else
    TempNA=Ini_NAx;
    Ini_NAx=-Ini_NAy;
    Ini_NAy=TempNA;
 end

% correct frequency coord of LED 
Ini_PixelShiftx=zeros(1,Length_MN);
Ini_PixelShifty=zeros(1,Length_MN);

for i=1:Length_MN
    pic_pos=i;
    
    Ini_PixelShiftx(pic_pos)=round(Ini_NAx(pic_pos)*Pixelsize*Nx);
    Ini_PixelShifty(pic_pos)=round(Ini_NAy(pic_pos)*Pixelsize*Ny);
end

%% Generate Complex or pure phase 3D object, and correspoding intensity bright field iamges

disp('Calculating slice-wise transfer functions...');
PTF_4D=single(zeros(Nx,Ny,length(Depth_Set),Length_MN));
ATF_4D=single(zeros(Nx,Ny,length(Depth_Set),Length_MN));

PTF_3D=single(zeros(Nx,Ny,length(Depth_Set)));
ATF_3D=single(zeros(Nx,Ny,length(Depth_Set)));


source=zeros(Nx,Ny);
figure
for i=1:Length_MN
    pic_pos=i
     
    source(Ini_PixelShiftx(pic_pos)+Nx/2+1,Ini_PixelShifty(pic_pos)+Ny/2+1)=1;
    for j=1:length(Depth_Set)
        deltaz=Depth_Set(j);

        [PTF,ATF]=model_TF(source,deltaz,Pixelsize,lambda,NA);
        PTF_3D(:,:,j)=PTF;
        ATF_3D(:,:,j)=ATF;  
    end
    PTF_4D(:,:,:,pic_pos)=PTF_3D;
    ATF_4D(:,:,:,pic_pos)=ATF_3D;

    source=zeros(Nx,Ny);        
end


%% Calculate Eq.(7) and Eq.(8) in paper
sum_PTF = 0;
sum_ATF = 0;

conj_PTF_Iten=0;
conj_ATF_Iten=0;

conj_term1=0;
conj_term2=0;

tic
for i=1:Length_MN
    pic_pos=i;

    Itmp = I_Raw(:,:,pic_pos);        
    Ihat_tmp = fft2((Itmp));

    sum_PTF=sum_PTF+abs(PTF_4D(:,:,:,pic_pos)).^2;
    sum_ATF=sum_ATF+abs(ATF_4D(:,:,:,pic_pos)).^2;

    conj_PTF_Iten=conj_PTF_Iten+conj(PTF_4D(:,:,:,pic_pos)).*repmat(Ihat_tmp,1,1,length(Depth_Set));
    conj_ATF_Iten=conj_ATF_Iten+conj(ATF_4D(:,:,:,pic_pos)).*repmat(Ihat_tmp,1,1,length(Depth_Set));

    conj_term1=conj_term1+conj(PTF_4D(:,:,:,pic_pos)).*ATF_4D(:,:,:,pic_pos);
    conj_term2=conj_term2+conj(ATF_4D(:,:,:,pic_pos)).*PTF_4D(:,:,:,pic_pos);
end

%% Repeat IDT algorithm

Normalized_term=(sum_PTF+Alpha).*(sum_ATF+Beta)-conj_term1.*conj_term1;
V_re = ((sum_ATF+Beta).* conj_PTF_Iten - conj_term1.* conj_ATF_Iten) ./ Normalized_term;
V_im = ((sum_PTF+Alpha).* conj_ATF_Iten - conj_term2.* conj_PTF_Iten) ./ Normalized_term;

toc
disp('Performing IDT spends...');

v_im = real(ifft2((V_im)));
v_re = real(ifft2((V_re)));
VV=-v_re+1i.*v_im;
RI = sqrt(VV./(k.^2) + n_Medium.^2);


