function [ xOut, A, weights ] = remOut( xIn, xOut, weights, method, alpha, scale, tol )
%Recursive function to remove outliers

%Calculate model
if strcmp(method,'rigidScale')
    %Rotation & translation
    [regParams, xEval, ErrorStats]=absor(xIn, xOut,'doTrans',true,'doScale',true,'weights',weights);

elseif strcmp(method,'rigid')
    %Rotation, translation, & scale
    [regParams, xEval, ErrorStats]=absor(xIn, xOut,'doTrans',true,'doScale',false,'weights',weights);

else
    %Affine
    inI=weights>tol; %Only include ones that have not been identified as outliers before
    A=padarray(xOut(:,inI),[1 0],1,'post')/padarray(xIn(:,inI),[1 0],1,'post'); %Fit to the data points
    xEval=A*padarray(xIn,[1 0],1,'post'); %Evaluate
    xEval(end,:)=[];
end

%Compare fit (xEval) to values (xOut)
fitDist=sqrt(sum((xEval-xOut).^2)); %Euclidean distance between points
outInd=false(size(fitDist));
outInd(weights~=0)=fitDist(weights~=0)>alpha.*std(fitDist(weights~=0));


if any(outInd(weights>tol))
    %If there are outliers where we're still counting the outlier a good deal, continue
    weights(outInd)=weights(outInd).*scale; %Downgrade the weight of outliers on the fit
    
    [ xOut, A, weights ] = remOut( xIn, xOut, weights, method,alpha, scale,tol ); %Call function again
else
    %If no more outliers, stop recursion
    xOut=xEval;
    
    if strcmp(method,'rigid') || strcmp(method,'rigidScale')
        A=regParams.M;
    end
    
end

