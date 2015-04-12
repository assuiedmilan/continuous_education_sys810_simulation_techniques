function [ oX, oY, oT ] = RK3(iCr, iBrs, iX0, iInput, iSystem, iNmax, iSampleTime)

wMatlabIndexBias = 1;

wA = iSystem{1};
wB = iSystem{2};
wC = iSystem{3};
wD = iSystem{4};

%Adding first rank of zeros for Brs coefficients if not already specified
%by user
if (all(iBrs(1,:)==0))
    wBrs = iBrs;
else
    wBrs = [zeros(1,size(iBrs,1));iBrs];
end

%Computing A(r)
wAr = sum(wBrs,2);

%Initialize
wRKNumberOfIterations = length(iCr);

wXp = zeros(wRKNumberOfIterations);
wUp = zeros(1,wRKNumberOfIterations);
oX = zeros(3,iNmax);
oX(:,1) = iX0;

for n = 0:iNmax-2
    
    wK = zeros(wRKNumberOfIterations);
    
    for wRank=1:wRKNumberOfIterations
        
        wSumBrsPerKs = zeros(wRKNumberOfIterations,1);
        
        for h=1:wRank-1
            wSumBrsPerKs(:) = wSumBrsPerKs(:) + wBrs(wRank,h).*wK(:,h);
        end
        
        wXp(:,wRank) = oX(:,n+wMatlabIndexBias) + iSampleTime*wAr(wRank)*(wSumBrsPerKs);
        wUp(:,wRank) = iInput(n+wMatlabIndexBias);
        wK(:,wRank)  = wA*wXp(:,wRank) + wB*wUp(:,wRank);
        
    end
    
    wSumCrPerKr = zeros(wRKNumberOfIterations,1);
    for h=1:wRKNumberOfIterations
        wSumCrPerKr(h) = wK(h,:)*iCr(:);
    end
    
    oX(:,n+wMatlabIndexBias+1) = oX(:,n+wMatlabIndexBias) + iSampleTime*(wSumCrPerKr);
    
end

oY = wC*oX + wD*iInput;

oT=linspace(0,iNmax*iSampleTime,iNmax);

end

