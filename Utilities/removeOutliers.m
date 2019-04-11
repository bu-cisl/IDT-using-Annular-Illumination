function [ freqDTh2, A ] = removeOutliers(freqDTh, freqDTh2, XYmid, method, alpha, scale, tol, DFI)
% Removes outliers by fitting a transformation from the original
% illumination angles (freqDTh) to the new, calibrated angles (freqDTh2);
% according to the different methods
% Outliers are ones that are far from the model fit; these are excluded
% from the model, then are replaced by modeled values

weights=ones(1,size(freqDTh,1)); %Weight all equally to start

if exist('DFI','var') && ~isempty(DFI)
    %If DFI is a nonempty variable, this means we want to remove outliers
    %for when DFI == 0 only (in the brightfield) and then extend the fit to
    %darkfield
    weights(logical(DFI))=0; %Do not use the darkfield to create fits
    
end

freqXY=pol2Cart(freqDTh, XYmid);
freqXY2=pol2Cart(freqDTh2,XYmid);

if strcmp(method,'affine') || strcmp(method,'projective')
    [transObj,ip1,~]=estimateGeometricTransform(freqXY(~DFI,:),freqXY2(~DFI,:),method);
    freqXY3=padarray(freqXY,[0 1],1,'post')*transObj.T; %Transform all values
    freqDTh3=cart2Pol(freqXY3(:,1:2), XYmid);%Convert back to polar
    outInd=~ismember(freqXY,ip1,'rows'); %Find the outliers + darkfield
else
    xIn=freqXY';
    xOut=freqXY2';
    [ xOut, A, weights ] = remOut( xIn, xOut, weights, method, alpha, scale, tol ); %Recursive function
    %When has been identified 2+ times as an outlier, replace the value with
    %the modeled value
    freqDTh3=cart2Pol(xOut', XYmid);
    outInd=weights'<tol; 
end
    
freqDTh2(outInd,:)=freqDTh3(outInd,:); %Only replace the outliers and the darkfield
  

end