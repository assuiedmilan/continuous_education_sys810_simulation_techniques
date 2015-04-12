function [ oX, oY, oT ] = RK3(iCr, iBrs, iX0, iInput, iSystem, iNmax, iSampleTime)

wMatlabIndexBias = 1;

wA = iSystem{1};
wB = iSystem{2};
wC = iSystem{3};
wD = iSystem{4};

oX = zeros(3,1);
oX(:,1) = iX0;

%Adding first rank of zeros for Brs coefficients if not already specified
%by user
if (all(iBrs(1,:)==0))
    wBrs = iBrs;
else
    wBrs = [zeros(1,size(iBrs,1));iBrs];
end

wRKNumberOfIterations = length(iCr);
for n = 0:iNmax-2
    
    wK = zeros(wRKNumberOfIterations);
    
    for wRank=1:wRKNumberOfIterations
        
        wSumK = zeros(wRKNumberOfIterations,1);
        
        for h=1:wRank-1
            wSumK(:) = wSumK(:) + wBrs(wRank,h).*wK(:,h)
        end
        
        wXp(:,wRank) = oX(:,n+wMatlabIndexBias) + iSampleTime*(wSumK);
        wUp(:,wRank) = iInput(n+wMatlabIndexBias);
        wK(:,wRank)  = wA*wXp(:,wRank) + wB*wUp(:,wRank);
        
    end
    
    wSumK = zeros(wRKNumberOfIterations,1);
    for h=1:wRKNumberOfIterations
        wSumK(h) = wK(h,:)*iCr(:);
    end
    
    oX(:,n+wMatlabIndexBias+1) = oX(:,n+wMatlabIndexBias) + iSampleTime*(wSumK);
    
end

oY = wC*oX + wD*iInput;

oT=linspace(0,iNmax*iSampleTime,iNmax);

end

