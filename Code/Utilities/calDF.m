function [ DFI ] = calDF( FI, XYmid, maxBk)
%FI, xI, yI, XYmid, radP, numImg)
%Identifies the darkfield images in the bunch - the calibration step won't
%work for these, so it's better to not try
%But then we can predict their positions based on the brightfield positions

%Get the DC
DC=squeeze(abs(FI(XYmid(2),XYmid(1),:)));
npix=numel(FI(:,:,1));

%Use the background values as well if they are available
if exist('maxBk','var') && ~isempty(maxBk)
    if size(maxBk,2)>size(maxBk,1)
        maxBk=maxBk';
    end
    X=[DC, maxBk];
else
    X=DC;
end

minBF=2000; %Heuristic value
if min(DC/npix)> minBF
    DFI=zeros(size(DC)); %Make sure don't cluster if all are brightfield
else
    Kidx=kmeans(X,2); %Cluster using k-means
    DFI=false(size(Kidx));

    %Assign darkfield to 1, brightfield to 0
    mDC=zeros(1,2);
    mDC(1)=mean(DC(Kidx==1)); %Mean of first cluster's DC
    mDC(2)=mean(DC(Kidx==2)); %Mean of second cluster's DC

    [~,mI]=min(mDC); %Get index of lower of two values

    DFI(Kidx==mI)=1; %Whichever cluster has lower mean DC is the darkfield
end

end

