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
        function oHandle = mDrawTimeseriesPlot(iThis,iValues,iTitle,iXlabel,iYlabel,varargin)
            
            if(isempty(varargin))
                
                oHandle = iThis.mProcessTimeseriesPlot(@plot,iValues,iTitle,iXlabel,iYlabel);
                
            else
                
                oHandle = iThis.mProcessTimeseriesPlot(str2func(varargin{1}),iValues,iTitle,iXlabel,iYlabel);
                
            end
            
        end
        
        function oHandle = mDrawStandardPlot(iThis,iValues,iPlot,varargin)           
               
                oHandle = iThis.mProcessStandardPlot(str2func(iPlot),iValues,varargin{:});
                
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
        
        function oHandle = mProcessTimeseriesPlot(iThis,iFuncHandle,iValues,iTitle,iXlabel,iYlabel)
            
            oHandle=figure();
            
            for i=1:length(iValues)
                
                iThis.mSetHold();
                
                iFuncHandle(iValues(i).Time,iValues(i).Data);
                wLegend = legend(get(legend(gca),'String'),iValues(i).Name);
                
            end
            
            set(wLegend, 'Interpreter', 'none');
            
            iThis.mProcessLabels(oHandle,iTitle,iXlabel,iYlabel);
            grid minor;
            iThis.mSaveDraw(oHandle);
        end
        
        function oHandle = mProcessStandardPlot(iThis,iFuncHandle,iValues,varargin)
            
            oHandle=figure();
            oHandle.Name = func2str(iFuncHandle);
            
            for i=1:2:length(iValues)
                
                iThis.mSetHold();
                
                if (i<length(iValues))
                    
                    iFuncHandle(iValues{i},iValues{i+1});
                else
                    iFuncHandle(iValues{i});
                end
                
            end
            
            iThis.mProcessLabels(oHandle,varargin{:});
            
            grid minor;
            iThis.mSaveDraw(oHandle);
        end
        
        function mProcessLabels(iThis,iHandle,varargin) %#ok<INUSL>
            
            set(0,'currentfigure',iHandle);
            
            if (length(varargin) == 1)
                title(varargin{1});
                iHandle.Name = varargin{1};
            elseif (length(varargin) == 2)
                xlabel(varargin{1});
                ylabel(varargin{2});
            elseif (length(varargin) == 3)
                title(varargin{1});
                iHandle.Name = varargin{1};
                xlabel(varargin{2});
                ylabel(varargin{3});
            end
        end
        
    end %Private methods
end %Class

