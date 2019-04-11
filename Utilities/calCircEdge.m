function [ freqDTh2, rad, pixDMat, pixD2Mat ] = calCircEdge (FIdivG, I, radP, freqDTh, XYmid, xI, yI, sigmaG, rScan, thScan, dScan, calRad, con, lambda)

%Row #1 = radius estimation, row #2 = center estimation
%rScan=[5 0.5]; %[+/- to scan, delta(r)]
%thScan=[10 1; 5 0.25];
%dScan=[20 1; 10 0.5];

%% Calibrate radius
centD2=freqDTh(:,1);
theta2=freqDTh(:,2);

imSz=size(FIdivG);
numImg=imSz(3);
imSz(3)=[];

if calRad
%Initialize
chRad=(rScan(1):-rScan(2):-rScan(1))'; %Range to scan radius in pixels
chTh=-thScan(1,1):thScan(1,2):thScan(1,1); %Angle of circle center from origin in degrees
chDi=-(dScan(1,1)/250)*imSz(1):dScan(1,2):(dScan(1,1)/250)*imSz(1); %Distance of circle center from origin
%Values based on a 20 pixel range on a 250x250 pixel image with circle
%radius = 60 pixels, so you (usually) don't have to change the parameter
%with changing image size

%% Determine radius

%Determine radius initialization
done=false;
checkDone=false;
numIter=0;

rad=radP;
radList=rad;
%medIList=0;

r_tol=[rScan(2)*1.1 rScan(2).*0.9];

PDindFR=1:(length(chRad)-1); %Allow any radius

fprintf('Finding radius: ')
while ~done
    %Keep track of iterations
    numIter=numIter+1;  
    
    %Print current image number
    fprintf('\r Iteration: %i; Radius: %f pixels',numIter,rad)
        
    radV=chRad + rad; %Radii to test
    
    %Select 20 images at random
    if numImg>20
        %Make sure none are the same
        imgI=checkSame(round(rand(1,20).*(numImg-1)) + 1,numImg);
    else
        imgI=1:numImg;
    end

    PDIfr=zeros(1,length(imgI));

    for jj=1:length(imgI)

        %Get the radial average of each selected image across the specified
        %radii
        centDV=chDi+centD2(imgI(jj)); %Make sure not going across the center (DC) of the image
        centDV(centDV<=0)=[];
        
        pixMean=radialAverage(FIdivG(:,:,imgI(jj)),rad,radV,chTh+theta2(imgI(jj)),centDV);
        %Differentiate across radius
        pixD=abs(diff(pixMean)); 

        %Select which radii to use - in this case, the one that has the max
        %signal
        [~,mPDI]=max(max(pixD(PDindFR,:),[],2));
        PDIfr(jj)=PDindFR(mPDI); %Record the index of the selected radius

    end
    
    medI=round(median(PDIfr)); %Take the median entry across the 20 images
    
    rad=radV(medI); %Select that radius
    
    %Store variables
    %medIList=[medIList medI];
    radList=[radList rad];
     
    %If not saying we're done
    if ~checkDone
        
        if abs(radList(end) - radList(end-1)) <= r_tol(1)
            %If change in radius is small enough, we're done
            checkDone=true;
            
        end
  
    else
        %If saying we are done
        radF= mean(radList(end-2:end)); %Average so not just taking random 20's consensus
        
        if abs(radF-radList(end)) <= r_tol(2)
            %If within tolerance
            rad=radF;
            done=true;
            
        else
            
            %If not in acceptable range or we have too much deviation
            %between iterations
            %(i.e. max should be only a few spaces off from radius we are positing)
            checkDone=false; %Otherwise, continue
            rad=radV(medI); %Not done, update the radius
            %Repeat with the new radius

        end
    end
  
    %Stop us from iterating endlessly
    if (numIter>=20 && ~checkDone) || rad <= 0
        rad=radP;
        done=true;
        %End it, put back to original predicted value
    end    
end
fprintf('\n');
NA_cal=(rad./con).*lambda;

fprintf('Radius = %.4f pixels; calibrated NA = %.4f\n',rad,NA_cal)

else
    rad=radP;
end
%% Circular edge detection

dR=0.5; 
%Initialize
chRad=(max(sigmaG+2*dR,3*dR):-dR:-dR)'; %Range to scan radius in pixels
chTh=-thScan(2,1):thScan(2,2):thScan(2,1); %Angle of circle center from origin in degrees
% chTh=permute(chTh,[3 1 2]); %Make theta stretch into 3rd dim
chDi=-(dScan(2,1)/250)*imSz(1):dScan(2,2):(dScan(2,1)/250)*imSz(1); %Distance of circle center from origin
%Scan to image size so can have constant dScan inputs, regardless of image
%size

radV=chRad + rad; 

numTh=length(chTh);
numD=length(chDi);
numR=length(chRad);

dRad1=radV(1:end-1)+diff(radV)./2; %Radii assoc. with first derivative
dRad2=radV(2:end-1);

PDI=length(dRad1); %First derivative will always be the last value in this formulation
PD2ind=find(dRad2>rad+dR); %Use the second derivative that's larger than the radius

%Set up angles defining circle
angleV=wrapTo180(0:180/64:360); %(in degrees), row vector
numA=length(angleV);

%Pre expand some of the matrices that don't depend on image
rMat=repmat(radV,[1 numA numTh numD]);
aMat=repmat(angleV,[numR 1 numTh numD]);
rcos=rMat.*cosd(aMat);
rsin=rMat.*sind(aMat);

%Storage
pixDMat=zeros(numTh,numD,numR-1,numImg);
pixD2Mat=zeros(numTh,numD,numR-2,numImg);
iMat=zeros(numImg,2);

maxPD=false(numTh,numD,numImg);
maxPD2=maxPD;
% maxBoth=maxPD;
% maxBoth2=maxPD;

maxBothxy=cell(numImg,2);
adjD=zeros(numImg,1);
adjTh=adjD;
%numPts=zeros(numImg,3);

filtXY=zeros(numImg,2,1); %Best x,y coord for each dataset (+overall best) & each metric

disp('Circular Edge Filter:')
fprintf('Evaluating image: ')
n=0;

for ii=1:numImg
    adjBndCount=0;
    %Print current image number
    fprintf(repmat('\b',1,n));
    strV=num2str(ii);
    fprintf(strV);
    n=numel(strV);
    
    %Select circle around the center point
    thetaV=wrapTo180(chTh+theta2(ii));
    distV=chDi+centD2(ii);
    
    done=false;
    while ~done

        %Deal with crossing origin 
        validD=true(size(distV));
        if any(distV<0)
            validD=distV>=0; %That is, want distV>=0
        end

        [arcMat, numA2]=calArc(rad, distV, thetaV, angleV); %Get arc
        arcMat=repmat(arcMat,[numR 1 1 numD]);

        %Rad in first dim, angle in 2nd
        %Put theta in 3rd dim, distance in 4th

        %Calculate xC, yC along axial spokes to try
        xC=repmat(distV,[numTh 1]).*cosd(repmat(thetaV',[1 numD]))+XYmid(1);
        yC=repmat(distV,[numTh 1]).*sind(repmat(thetaV',[1 numD]))+XYmid(2);

        %Calculate the x,y coordinates of the circle edges
        xR=reshape(rcos(arcMat),[numR numA2 numTh numD])+repmat(permute(xC,[3 4 1 2]),[numR numA2 1 1]);
        yR=reshape(rsin(arcMat),[numR numA2 numTh numD])+repmat(permute(yC,[3 4 1 2]),[numR numA2 1 1]);

        %Interpolate the circle values at these locations
        pixVal=interp2(xI,yI,FIdivG(:,:,ii),xR,yR,'spline');
        %pixMat=reshape(pixVal,[numR numA numTh numD]);

        %Derivatives
        pixMean=squeeze(mean(pixVal,2)); %Mean value around rim

        %1st derivative
        pixD=diff(pixMean);

        %2nd derivative
        pixD2=diff(pixMean,2);

        %Select which radii to use for the second derivative
        %Choose the one with max value, but only larger than the actual radius
        %Only use radii larger than the radius
        [~,mPD2I]=max(max(pixD2(PD2ind,:),[],2));
        PD2I=PD2ind(mPD2I);

        %Save it
        iMat(ii,:)=[PDI,PD2I];

        %Store
        pixDMat(:,:,:,ii)=permute(pixD,[2 3 1]);
        pixD2Mat(:,:,:,ii)=permute(pixD2,[2 3 1]);

        %1st Derivative
        pd=squeeze(pixD(PDI,:,:)); %Select radius
        pd(:,~validD)=0;%Deal with case where distV<0 (so don't get circle on other side)
        mxL=max(pd(:));
        sL=std2(pd);
        maxPD(:,:,ii)=(pd>=mxL-0.1.*sL); %Select region around max

        %2nd Derivative
        pd2=squeeze(pixD2(PD2I,:,:)); %Select radius
        pd2(:,~validD)=0;
        mxL2=max(pd2(:));
        sL2=std2(pd2);
        maxPD2(:,:,ii)=(pd2>=mxL2-0.25.*sL2); %Select region around max


        %Combine them both into one map
        maxBoth=maxPD(:,:,ii) & maxPD2(:,:,ii);

        if ~any(any(maxBoth))
            %If there's no overlap between the two, just use the 1st derivative
            %map
            maxBoth=maxPD(:,:,ii);
        end
        maxBothxy(ii,:)=[{xC(maxBoth)},{yC(maxBoth)}]; %Save the possible coordinates

        %Check if the only triggered pixels are on the extreme edges of the
        %search area - means need to change search area
        [thI,dI]=find(maxBoth);
        unDi=unique(dI); unTh=unique(thI);

        adjBnd=false;
        if all(unDi<=2) || all(unDi>=numD-1)
            %If all distances are at the edge 
            if unDi(1)<=2
                adjD(ii)=distV(unDi(1)); %Center on one of them
            else
                adjD(ii)=distV(unDi(end)); %Center on one of them         
            end
            distV=adjD(ii)+chDi;
            adjBnd=true; %Adjust bounds is true --> iterate
        end

        if all(unTh<=2) || all(unTh>=numTh-1)
            %If all thetas are at the edge --> iterate
            if unTh(1)<=2
                adjTh(ii)=thetaV(unTh(1)); %Center on 
            else
                adjTh(ii)=thetaV(unTh(end)); %Center on one
            end
            thetaV=wrapTo180(chTh+adjTh(ii));
            adjBnd=true;
        end

        if adjBndCount >= 20
            %Stop from endlessly iterating
            adjBnd=false;
        end

        if ~adjBnd
            %If not adjusting bounds, conclude
            done=true;    
            tempPt=[maxBothxy{ii,1},maxBothxy{ii,2}]; %Extract x,y coordinates

            [tempErr]=imageErr(I(:,:,ii),tempPt(:,1),tempPt(:,2),rad,xI,yI);
        %tempErr=-imageErr(I(:,:,ii),tempPt(:,1),tempPt(:,2),rad,xI,yI);
            %Select the circle center that gives minimum RMSE error
            [~,mI]=min(tempErr,[],1);
            filtXY(ii,1:2)=tempPt(mI,1:2);
        else
            %Iterate if adjusting bounds
            fprintf('\nAdjusting image %i to distance %f pixels, theta %f degrees.\n',ii,adjD(ii),adjTh(ii))
            fprintf('Evaluating image: ');
            n=0;
            adjBndCount=adjBndCount+1;
        end
    end
    

end

freqDTh2=cart2Pol(filtXY,XYmid);

fprintf('\n')

end
