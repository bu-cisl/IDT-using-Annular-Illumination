function [ XYopp ] = convOpp( freqXY, XYmid )
%Get the opposite circle centers

XYopp=2.*repmat(XYmid,[size(freqXY,1) 1])-freqXY;

end

