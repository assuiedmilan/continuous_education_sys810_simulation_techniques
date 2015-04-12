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
    
    wSumK = zeros(3,1);
    for h=1:1
        wSumK(:,h) = iBrs(1,h).*wK(:,h);
    end
    
    wXp(:,1) = oX(:,n+wMatlabIndexBias) + iSampleTime*(sum(wSumK,2));
    wUp(:,1) = iInput(n+wMatlabIndexBias);
    
    wK(:,2)  = wA*wXp(:,1) + wB*wUp(:,1);
    
    wSumK = zeros(3,1);
    for h=1:2
        wSumK(:,h) = iBrs(2,h).*wK(:,h);
    end
    
    wXp(:,2) = oX(:,n+wMatlabIndexBias) + iSampleTime*(sum(wSumK,2));
    wUp(:,2) = iInput(n+wMatlabIndexBias);
    
    wK(:,3) = wA*wXp(:,2) + wB*wUp(:,2);
    
    wSumK = zeros(3,1);
    for h=1:3
        wSumK(h) = wK(h,:)*iCr(:);
    end
    
    oX(:,n+wMatlabIndexBias+1) = oX(:,n+wMatlabIndexBias) + iSampleTime*(wSumK);
    
end

oY = wC*oX + wD*iInput;

oT=linspace(0,iNmax*iSampleTime,iNmax);

end

