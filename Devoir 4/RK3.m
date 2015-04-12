function [ oX, oY, oT ] = RK3(iCr, iBrs, iX0, iInput, iSystem, iNmax, iSampleTime)

wMatlabIndexBias = 1;

wA = iSystem{1};
wB = iSystem{2};
wC = iSystem{3};
wD = iSystem{4};

oX = zeros(3,1);
oX(:,1) = iX0;

for n = 0:iNmax-2
    
    wK(:,1)  = wA*oX(:,n+wMatlabIndexBias) + wB*iInput(n+wMatlabIndexBias);
        
    wXp(:,1) = oX(:,n+wMatlabIndexBias) + iSampleTime*(sum(iBrs(1,1:1).*wK(1:1)));    
    wUp(:,1) = iInput(n+wMatlabIndexBias);
    
    wK(:,2)  = wA*wXp(:,1) + wB*wUp(:,1);

    wXp(:,2) = oX(:,n+wMatlabIndexBias) + iSampleTime*(sum(iBrs(2,1:2).*wK(1:2)));
    wUp(:,2) = iInput(n+wMatlabIndexBias);
    
    wK(:,3) = wA*wXp(:,2) + wB*wUp(:,2);
    
    wKsum = zeros(size(wK,1),1);
    for h=1:size(wK,1)
        wKsum(h,1) = sum(iCr(:).*wK(:,h));
    end
    
    oX(:,n+wMatlabIndexBias+1) = oX(:,n+wMatlabIndexBias) + 1/sum(iCr)*iSampleTime.*wKsum;
    
end

oY = wC*oX + wD*iInput;

oT=linspace(0,iNmax*iSampleTime,iNmax);

end

