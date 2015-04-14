classdef RK3 < handle
    
    methods(Static)
        
        function [ oX,  oT ] = sProcessRungeKutta(iCr, iBr, iX0, iInput, iSystem, iNmax, iSampleTime, iDXhandle)
            
            wMatlabIndexBias = 1;
            
            %Adding first rank of zeros for Brs coefficients if not already specified
            %by user
            if (all(iBr(1,:)==0))
                wBrs = iBr;
            else
                wBrs = [zeros(1,size(iBr,1));iBr];
            end
            
            %Computing A(r)
            wAr = sum(wBrs,2);
            
            %Initialize
            wRKNumberOfIterations = length(iCr);
            
            oX = zeros(3,iNmax);
            oT = zeros(1,iNmax);
            oX(:,1) = iX0;
            oT(1) = 0;
            
            for n = 0:iNmax-2
                
                wK = zeros(wRKNumberOfIterations);
                
                for wRank=1:wRKNumberOfIterations
                    
                    wSumBrsPerKs = zeros(wRKNumberOfIterations,1);
                    
                    for h=1:wRank-1
                        wSumBrsPerKs(:) = wSumBrsPerKs(:) + wBrs(wRank,h).*wK(:,h);
                    end
                    
                    wK(:,wRank)  = iDXhandle(iSampleTime*(n + wAr), oX(:,n+wMatlabIndexBias) + iSampleTime*(wSumBrsPerKs), iInput(n+wMatlabIndexBias));
                    
                end
                
                wSumCrPerKr = zeros(wRKNumberOfIterations,1);
                for h=1:wRKNumberOfIterations
                    wSumCrPerKr(h) = wK(h,:)*iCr(:);
                end
                
                oX(:,n+wMatlabIndexBias+1) = oX(:,n+wMatlabIndexBias) + iSampleTime*(wSumCrPerKr);
                oT(n+wMatlabIndexBias+1) = (n+1)*iSampleTime;
                
            end
            
        end
        
        function oCr = sComputeCrOrder3(iCr3,iBr,varargin)
            
            if(isempty(varargin))
                %Precision 3 Kutta's method
                g = 6;
            else
                %Precision 2 Kutta's method
                g = varargin{1};
            end
            
            %Kutta's third-order precision 2 method: g=6 => precision 3 Kutta's method
            %0  |  0         0                        0
            %a2 |  b21       0                        0
            %a3 |  b31       b32                      0
            %------------------------------------------------------
            %    | [1-c2-c3]  [1/a2 * (1/2 - c3*a2)] [1/g *c3]
            
            wAr = sum(iBr,2);
            
            wC3  = 6/g*iCr3;
            wC2  = 1/wAr(1) * (1/2 - wC3*wAr(2));
            wC1  = 1 - wC2 - wC3;
            
            oCr = [wC1,wC2,wC3];
            
        end
    end
    
end
