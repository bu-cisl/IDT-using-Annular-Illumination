function [ freqUV ] = fitToRigidTransform( freqUV, freqUV_design, method, alpha, scale, tol )
%Fits to rigid transform in order to remove outliers
%Does NOT completely wipe out the structure of the data - just removes the
%complete outliers and finds the best fit to the remaining points. That fit
%is then applied to the outliers ONLY, not to all data points

weights=ones(1,size(freqUV,1)); % Weight all equally

if strcmp(method,'affine') || strcmp(method,'projective')
    [transObj,ip1,~] = estimateGeometricTransform(freqUV_design, freqUV, method);
    freqUV3=padarray(freqUV_design,[0 1],1,'post')*transObj.T; %Transform all values
    outInd=~ismember(freqUV_design,ip1,'rows'); %Find the outliers + darkfield
else
    xIn=freqUV_design';
    xOut=freqUV';

    [ xOut, A, weights ] = remOut( xIn, xOut, weights, method, alpha, scale, tol ); %figH, figH2
    %When has been identified 2+ times as an outlier, replace the value with
    %the modeled value
    freqUV3=xOut';
    outInd=weights'<tol; 
end
    
freqUV(outInd,:)=freqUV3(outInd,:); %Only replace the outliers
end

