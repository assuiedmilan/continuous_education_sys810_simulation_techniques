function [ oX, oY, oT ] = RK3(iCr, iBrs, iX0, iInput, iSystem, iNmax, iSampleTime)

wMatlabIndexBias = 1;


wA = iSystem{1};
wB = iSystem{2};
wC = iSystem{3};
wD = iSystem{4};

oX = zeros(3,1);
oX(:,1) = X0;

%Compute Ar:
wAr = sum(iBrs,2);

for n = 0:iNmax-1
    
    wK(1)  = wA*oX(:,n+wMatlabIndexBias) + wB*iInput(n+wMatlabIndexBias);
        
    wXp(1) = oX(:,n+wMatlabIndexBias) + iSampleTime(sum(iBrs(2,1:1).*wK(1:1)));    
    wUp(1) = iInput(n+wMatlabIndexBias);
    
    wK(2)  = wA*wXp(1) + wB*wUp(1);

    wXp(2) = oX(:,n+wMatlabIndexBias) + iSampleTime(sum(iBrs(3,1:2).*wK(1:2)));
    wUp(2) = iInput(n+wMatlabIndexBias);
    
    wK(3) = wA*wXp(2) + wB*wUp(2);
    
    oX(:,n+wMatlabIndexBias+1) = oX(:,n+wMatlabIndexBias) + 1/sum(wAr)*iSampleTime*(sum(wAr.*wK));
    
end

oY = wC*oX + wD*iInput;

oT=linspace(0,iSimulationTime,iSimulationTime*(1/iSampleTime));

end

