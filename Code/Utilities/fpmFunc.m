function [objR, Pup, err, scale, freqUV, time] = ptychReconstruct_crop_v2( I, nObj, freqUV, opts)
%AlterMin Implements alternating minimization sequentially on a stack of
%measurement I (n1 x n2 x nz). It consists of 2 loops. The main loop updates
%the reconstruction results O and P. The inner loop applies projectors/minimizers
%P1 and P2 on each image I and steps through the entire dataset.


%% Default options
if nargin < 4
    % default values
    opts.tol = 1;
    opts.maxIter = 50;
    opts.minIter = 3;
    opts.monotone = 1;
    opts.display = 0;
    opts.saveIterResult = 0;
    opts.out_dir = [];
    opts.O0 = Ft(sqrt(I(:,:,1)))/r0;
    opts.O0 = padarray(opts.O0,(nObj-NsampR)/2);
    opts.P0 = ones(NsampR);
    opts.OP_alpha = 1;
    opts.OP_beta = 1;
    opts.mode = 'real';
    opts.scale = ones(Nled,1);
    opts.H0 = ones(NsampR);
    opts.position_calibration_type = 0;
    opts.calbratetol = 1e-1;
    opts.F = @(x) ifftshift(fft2(fftshift(x)));
    opts.Ft = @(x) ifftshift(ifft2(fftshift(x)));
    opts.scanRange=1;
    opts.radialPenalty=0;
    opts.gradIter=1;
    opts.freqUV_design=freqUV;
    opts.fitRigidTransform=false;
    opts.transMethod='rigidScale';
    opts.transAlpha=2;
    opts.transScale=0.5;
    opts.transTol=0.05;
end

%% Operators & derived constants
% size of measurement
NsampR = size(I); %Size of real space
Nimg = NsampR(3);
NsampR(3) = [];

cen0 = round((nObj+1)/2); %Coordinate of center of Fourier space

%% Operators
row = @(x) x(:).';

% operator to crop region of O from proper location at the O plane
downsamp = @(x,cen) x(cen(2)-floor(NsampR(1)/2):cen(2)-floor(NsampR(1)/2)+NsampR(1)-1,...
    cen(1)-floor(NsampR(2)/2):cen(1)-floor(NsampR(2)/2)+NsampR(2)-1);
%cen is in (x,y)
%For the given image x, crop out a rectangle of NsampR(1) x NsampR(2) row, col
%extent around center cen

%% Initialization
%Copy some parameters
Ps = opts.Ps;
H0 = opts.H0; 
F = opts.F;
Ft = opts.Ft;

Pup = opts.P0; opts.P0 = 0; %Pupil function in Fourier domain
objR = opts.O0; opts.O0 = 0; %Object in real space
ObjF=F(objR); %Init in Fourier space
%ObjF=padarray(F(objR),[(nObj(1)-NsampR(1))/2,(nObj(2)-NsampR(2))/2]);
%ps=opts.dpix_mF; %Pixel size in recon image
con=opts.con;

err1 = inf;
err2 = 50;
err = zeros(1,opts.maxIter);
iter = 0;
scale = opts.scale;

%Spatial coordinates for object space
x_obj = opts.x_obj;
y_obj = opts.y_obj;

%Debugging parameters
% ObjF_init=zeros([size(x_obj),Nimg]);
% dPsi=zeros([size(Pup),Nimg]);
% Psi=dPsi;
% ObjF_save=ObjF_init;
% objR_save=ObjF_init;
% objRsave=zeros([size(opts.O0) opts.maxIter+1]);

%% Display init
fprintf('|  Iteration  |   Image   |      Cost     |  LED Pos. Step |\n');
fprintf('+-------------+-----------+---------------+----------------+\n');
fprintf('|     %3d     |    -/-    |       -       |     %.3f      |\n',iter,0);

if opts.display
    if strcmp(opts.mode,'real')
        o = objR;
    elseif strcmp(opts.mode,'fourier')
        o = F(objR);
    end
    f1 = figure(88);
    if ~contains(opts.position_calibration_type, 'online')
        subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
        title('Object Amp.');
        subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
        title('Object Phase');
        subplot(223); imagesc(abs(Pup)); axis image; colormap gray; colorbar;
        title('Pupil Amp');
        subplot(224); imagesc(angle(Pup).*abs(Pup)); axis image; colorbar;
        title('Pupil Phase');
    else
        subplot(2,4,1); imagesc(abs(o)); axis image; colormap gray; colorbar;
        title('Object Amp.');
        subplot(2,4,2); imagesc(angle(o)); axis image; colormap gray; colorbar;
        title('Object Phase');
        subplot(2,4,5); imagesc(abs(Pup)); axis image; colormap gray; colorbar;
        title('Pupil Amp');
        subplot(2,4,6); imagesc(angle(Pup).*abs(Pup)); axis image; colorbar;
        title('Pupil Phase');

        subplot(2,4,[3,4,7,8]); scatter(freqUV(:,1), freqUV(:,2)); axis image
        title('LED Positions');
        xlabel('NA_x'); ylabel('NA_y');
        xlim([-1,1]); ylim([-1,1]);
    end
    drawnow;
end

if opts.saveIterResult
    export_fig(f1,[opts.out_dir,'\R_',num2str(iter),'.png'],'-m4');
end


%% Main 
% stopping criteria: when relative change in error falls below some value,
% can change this value to speed up the process by using a larger value but
% will trading off the reconstruction accuracy
% error is defined by the difference b/w the measurement and the estimated
% images in each iteration

% Convert from 1/um center to corresponding crop region
% **** Check this conversion
Pupilshifty = round(freqUV(:,2)*con);
Pupilshiftx = round(freqUV(:,1)*con);

%Corresponding crop regions
YXmid=floor(nObj./2)+1;
hfSz=floor(NsampR(1)/2);
cropR=zeros(Nimg,4);
cropR(:,1) = YXmid(1)+Pupilshifty-hfSz; %y start
cropR(:,2) = YXmid(1)+Pupilshifty+hfSz-1; %y end
cropR(:,3) = YXmid(2)+Pupilshiftx-hfSz; %x start
cropR(:,4) = YXmid(2)+Pupilshiftx+hfSz-1; %x end

%Simulated annealing parameters
freqXY=cropR(:,[3 1])+hfSz;
sp0 = max(abs(freqXY(1,:)-freqXY(2,:))); %Only search between adjacent illum

%Corresponding to freqUV=(0,0)
cropR_00 = [YXmid(1)-hfSz, ... %y start
        YXmid(1)+hfSz-1, ... %y end
        YXmid(2)-hfSz, ... %x start
        YXmid(2)+hfSz-1]; % end

Objcrop = zeros([NsampR Nimg]);  

tic
while abs(err1 - err2) > opts.tol && iter < opts.maxIter
    err1 = err2;
    err2 = 0;
    iter = iter+1; 
    %objRsave(:,:,iter)=objR;
    fprintf('Evaluating image: ')
    strV=blanks(4);
    fprintf('%s',strV) 
    for ii = 1:Nimg           
        %Print current image number
        str=num2str(ii);
        strV(1:length(str))=str;
        fprintf('\b\b\b\b%s',strV)
        
        % Crop measured intensity
        I_mea = I(:,:,ii);
        Objfcrop = ObjF(cropR(ii,1):cropR(ii,2),cropR(ii,3):cropR(ii,4));
        ObjfcropP = Objfcrop.*Pup;
        ObjcropP = Ft(ObjfcropP);
        Objfup = F(sqrt(I_mea).*ObjcropP./(abs(ObjcropP+eps)));
        
        % Pupil Update
        Pup = Pup + abs(Objfcrop) .* conj(Objfcrop) .* (Objfup - ObjfcropP) / max(abs(ObjF(:))) ./ (abs(Objfcrop).^2 + opts.OP_beta) .* Ps;
        
        % (Cropped) Object update
        ObjF(cropR(ii,1):cropR(ii,2),cropR(ii,3):cropR(ii,4)) = ObjF(cropR(ii,1):cropR(ii,2),cropR(ii,3):cropR(ii,4)) + ...
            abs(Pup) .* conj(Pup) .* (Objfup - ObjfcropP) / max(abs(Pup(:))) ./ (abs(Pup).^2 + opts.OP_alpha);
        
        % Update sequential 
        Objcrop(:,:,ii) = Ft(ObjF(cropR(ii,1):cropR(ii,2),cropR(ii,3):cropR(ii,4)).*Pup);

        % Position correction
        poscost = @(ss) sum(sum((sqrt(abs(Ft(downsamp(ObjF, round(ss)).*Pup.*H0)).^2)-sqrt(I_mea)).^2));
        
        if contains(opts.position_calibration_type, 'online')
            [dk_x, dk_y] = findLedNaError(I_mea, ObjF, Pup, cropR(ii,:), opts.scanRange, opts.radialPenalty, opts.gradIter);
            cropR(ii,1:2)=cropR(ii,1:2) + dk_y;
            cropR(ii,3:4)=cropR(ii,3:4) + dk_x;
        end
    end
    fprintf('\n')
    
    objR = Ft(ObjF);

    %% Compute error
    % record the error and can check the convergence later.
    err(iter)= sum(sum(sum(((I-abs(Objcrop).^2).^2))));
    err2=err(iter);
%     fprintf('| %2d   | %.2e |\n',iter,err2);

    
    if opts.saveIterResult
        export_fig(f1,[opts.out_dir,'\R_',num2str(iter),'.png'],'-m4');
        %     saveas(f2,[opts.out_dir,'\Ph_',num2str(iter),'.png']);
    end
    
    if opts.monotone&&iter>opts.minIter
        if err2>err1
            break;
        end
    end
    
    %Convert back to frequency coordinates
    Pupilshifty = cropR(:,1) + hfSz - YXmid(1);
    Pupilshiftx = cropR(:,3) + hfSz - YXmid(2);
    freqUV_prev = freqUV;
    freqUV=[Pupilshiftx Pupilshifty] ./ con; %Conversion from pixels to freq space
    
    % Perform rigid fit of calibrated LED positions if settings indicate
    % this
    if contains(opts.position_calibration_type, 'online') ...
            && contains(opts.position_calibration_type, 'rigid') ...
            && ~contains(opts.position_calibration_type, 'nonrigid')
        
        [ freqUV ] = fitToRigidTransform(freqUV, opts.freqUV_design, ...
                                         opts.transMethod, ...
                                         opts.transAlpha, ...
                                         opts.transScale, ...
                                         opts.transTol);
        
        % Determine pupil shift              
        Pupilshifty = round(freqUV(:,2) * con);
        Pupilshiftx = round(freqUV(:,1) * con);

        % Determine corresponding crop regions
        cropR(:,1) = YXmid(1) + Pupilshifty - hfSz; %y start
        cropR(:,2) = YXmid(1) + Pupilshifty + hfSz - 1; %y end
        cropR(:,3) = YXmid(2) + Pupilshiftx - hfSz; %x start
        cropR(:,4) = YXmid(2) + Pupilshiftx + hfSz - 1; %x end
    end
    
    % Display results
    if opts.display
        if strcmp(opts.mode,'real')
            o = objR;
        elseif strcmp(opts.mode,'fourier')
            o = F(objR);
        end
        f1 = figure(88);
        if ~contains(opts.position_calibration_type, 'online')
            % Display object and pupil only
            subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
            title('Object Amp.');
            subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
            title('Object Phase');
            subplot(223); imagesc(abs(Pup)); axis image; colormap gray; colorbar;
            title('Pupil Amp');
            subplot(224); imagesc(angle(Pup).*abs(Pup)); axis image; colorbar;
            title('Pupil Phase');
        else
            % Display object, pupil, and LED positions
            subplot(2,4,1); imagesc(abs(o)); axis image; colormap gray; colorbar;
            title('Object Amp.');
            subplot(2,4,2); imagesc(angle(o)); axis image; colormap gray; colorbar;
            title('Object Phase');
            subplot(2,4,5); imagesc(abs(Pup)); axis image; colormap gray; colorbar;
            title('Pupil Amp');
            subplot(2,4,6); imagesc(angle(Pup).*abs(Pup)); axis image; colorbar;
            title('Pupil Phase');

            subplot(2,4,[3,4,7,8]); scatter(freqUV(:,1), freqUV(:,2)); axis image
            xlabel('NA_x'); ylabel('NA_y');
            xlim([-1,1]); ylim([-1,1]);
            title('LED Positions');
        end
        drawnow;
    end
    
    fprintf('|     %3d     |    -/-    |    %.2e   |     %.3f      |\n',iter,err2,norm(freqUV_prev - freqUV));
end


end

