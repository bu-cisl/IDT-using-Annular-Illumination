function [ vecOut ] = checkSame( vecIn, maxVal )
%Check if any of the values are the same. If so, change one of them.

[vecx,vecy]=meshgrid(vecIn); %Check if two values are the same 
difVec=vecx==vecy;
difVec(logical(eye(size(difVec))))=0;
difVec=sum(difVec);
if any(difVec(:)) %If any are the same
    loc=find(difVec,1); %Find them
    vecIn(loc)=vecIn(loc)+1; %Add 1 to the matching one
    if vecIn(loc)>maxVal %If we've gone beyond bounds, reset
        vecIn(loc)=1;
    end
    vecOut=checkSame(vecIn, maxVal); %Recursion
else    
    vecOut=vecIn; %End the recursion
end

end

