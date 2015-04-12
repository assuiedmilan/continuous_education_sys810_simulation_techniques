classdef RK3 < handle
    
    methods(Static)
        
        function [ oX, oY, oT ] = sProcessRungeKutta(iCr, iBrs, iX0, iInput, iSystem, iNmax, iSampleTime)
            
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
        
        function oHandle = sAddButcherTable(iHandle,iCr, iBrs)
            
            set(0, 'currentfigure', iHandle);
            hold all;
            subplot(2,1,1);
            
            %Adding first rank of zeros for Brs coefficients if not already specified
            %by user
            if (all(iBrs(1,:)==0))
                wBrs = iBrs;
            else
                wBrs = [zeros(1,size(iBrs,1));iBrs];
            end
            %Add one last column of not significants zeros for above last
            %Cr            
            wBrs = [wBrs,zeros(size(wBrs,1),1)];
            
            %Adding one not significant 0 to Crs
            wCrs = [0,iCr];
            
            %Computing A(r)
            wAr = sum(wBrs,2);
            
            wData = [wAr,wBrs];
            wData = [wData;wCrs];
            
            wLineName = cell(length(wAr)+1,1);
            wLineName{length(wAr)+1} = 'Cr';

            uitable('Data', wData, 'ColumnName', {'Ar', 'Brs'}, 'RowName', wLineName, 'Position', [20 20 500 150]);
                       
            
            oHandle = iHandle;
            
                        
        end
        
    end
    
end
