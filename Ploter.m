classdef Ploter < handle
    
    properties (SetAccess = private, GetAccess = private)
        
        mPaperPos;
        mPaperSize
        mHoldAll;
        mSaveAfterDraw;
        
    end
    
    methods %Public
        
        %Constructors
        function oInstance = Ploter(iPaperPos, iPaperSize, varargin)
            
            oInstance.mPaperPos = iPaperPos;
            oInstance.mPaperSize = iPaperSize;
            
            oInstance.mHoldAll = true;
            oInstance.mSaveAfterDraw = true;
            
        end
        
        %Public methods
        function oHandle = mDrawStandardPlot(iThis,iTimeSeries,iXLabel,iYLabel,iTitle)
            
            oHandle=figure();
            oHandle.Name = iTitle;
            
            for i=1:length(iTimeSeries)
                
                iThis.mSetHold();
                
                plot(iTimeSeries(i).Time,iTimeSeries(i).Data);
                wLegend = legend(get(legend(gca),'String'),iTimeSeries(i).Name);
                
            end
            
            set(wLegend, 'Interpreter', 'none');
            xlabel(iXLabel);
            ylabel(iYLabel);
            title(iTitle);
            grid minor;
            iThis.mSaveDraw(oHandle);
        end
        
        %Accessors
        
        function mSetPaperPos(iThis,iPaperPos)
            iThis.mPaperPos = iPaperPos;
        end
        
        function mSetPaperSize(iThis,iPaperSize)
            iThis.mPaperSize = iPaperSize;
        end
        
    end %Public methods
    
    methods (Access = private)
        
        function mSetHold(iThis)
            
            if (iThis.mHoldAll)
                hold all;
            end
            
        end
        
        function mSaveDraw(iThis,iHandle)
            
            if (iThis.mSaveAfterDraw)
                iThis.mProcessSaveDraw(iHandle);
            end
            
        end
        
        function mProcessSaveDraw(iThis,iHandle)
            
            set(iHandle, 'PaperPosition', iThis.mPaperPos);
            set(iHandle, 'PaperSize', iThis.mPaperSize);
            saveas(iHandle, iHandle.Name, 'pdf');
            
        end
        
    end %Private methods
end %Class

