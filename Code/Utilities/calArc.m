function [ arcMat, numA2 ] = calArc( rad, distV, thetaV, angleV)
    %Select arc so don't use center region
    %Only let vary with theta so have the same number of angles
    numTh=length(thetaV);
    numA=length(angleV);
    
    if mean(distV)<=rad
        phi=acosd(mean(distV)./rad);
    else
        phi=0;
    end
    thBnd2=(thetaV)+phi-180; %Lower bound
    thBnd3=2.*(thetaV) - thBnd2; %Upper bound

    %Deal with angle wrapping issues
    th2w=wrapTo180(thBnd2);    
    th3w=wrapTo180(thBnd3);
    angleA=repmat(angleV,[numTh 1])>=repmat(th2w',[1 numA]);
    angleB=repmat(angleV,[numTh 1])<=repmat(th3w',[1 numA]);

    %Define logical index matrix to select angles
    arcMat=angleA&angleB;
    wrapLoc=repmat(thBnd2',[1 numA])<-180 |repmat(thBnd3',[1 numA])>180;
    arcMat(wrapLoc)=angleA(wrapLoc)|angleB(wrapLoc);
    
    %Make sure number of angles requested across angles is the same (data
    %size)
    numAcap=sum(arcMat,2);
    numA2=mode(numAcap);
    difA=numAcap-numA2; %Need to add or subtract angles
    
    if any(difA<0)
        %Need to add angles
        ind=find(difA<0);
        for kk=1:length(ind) %Loop across all vectors
            numAdd=abs(difA(ind(kk))); %Number to add
            temp=arcMat(ind(kk),:); 
            dT=diff(temp); %Find edges
            even=floor(numAdd/2); %Find even number to add
            odd=mod(numAdd,2); %Find if need to add additional 1
            
            up=find(dT<0); %Upper end (where go 1-->0)
            bot=find(dT>0);%Lower end (where go 0-->1)
            
            if isempty(up)
                up=length(temp); %If none, first in vector is 0 , last is 1 
            end
            if isempty(bot)
                bot=1;  %If none, first is 1, last is 0
            end
            
            upI=up+1:up+even+odd; %Add odd and half of the even to the upper end
            upI(upI>length(temp))=upI(upI>length(temp))-length(temp); %Wrap to front if go off the end
            
            botI=bot:-1:bot-even+1; %Add half of even to front
            botI(botI<1)=botI(botI<1)+length(temp); %Wrap around
            
            temp([upI botI])=1; %Use indices to change value
            arcMat(ind(kk),:)=temp; %Re-store
        end
    end
    if any(difA>0)
        %Need to delete angles
        ind=find(difA>0);
        for kk=1:length(ind)
            numDel=abs(difA(ind(kk)));%Loop across all vectors
            temp=arcMat(ind(kk),:);%Number to delete
            dT=diff(temp);%Find edges
            even=floor(numDel/2);
            odd=mod(numDel,2);
            
            up=find(dT<0);
            bot=find(dT>0);
                        
            if isempty(up)
                up=length(temp);  
            end
            if isempty(bot)
                bot=1;  
            end
            
            upI=up:-1:up-even-odd+1;
            upI(upI<1)=upI(upI<1)+length(temp);
            
            botI=bot+1:bot+even;
            botI(botI>length(temp))=botI(botI>length(temp))-length(temp);
            
            temp([upI botI])=0;
            arcMat(ind(kk),:)=temp;
        end
    end

    arcMat=permute(arcMat,[3 2 1]);

end

