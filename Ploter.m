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
        
        %Public drawers
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
        
        function mSetSaveAfterDraw(iThis,iValue)
            iThis.mSaveAfterDraw = iValue;
        end
        
        %Save methods
        function mProcessSaveDraw(iThis,iHandle)
            
            set(iHandle, 'PaperPosition', iThis.mPaperPos);
            set(iHandle, 'PaperSize', iThis.mPaperSize);
            %waitfor(iHandle);
            saveas(iHandle, iHandle.Name, 'pdf');
            
        end
        
        function mProcessSaveModel(iThis,iBlock)
            
            saveas(get_param(iBlock,'Handle'), strcat('SimBlock_',get_param(iBlock,'Name')), 'pdf');
            
        end
        
    end %Public methods
    
    methods (Access = private)
        
        function mSetHold(iThis)
            
            if (iThis.mHoldAll)
                hold all;
            end
            
        end

        function mSetSaveDraw(iThis,iValue)
            
            iThis.mSaveAfterDraw = iValue;
        end
        
        function mSaveDraw(iThis,iHandle)
            
            if (iThis.mSaveAfterDraw)
                iThis.mProcessSaveDraw(iHandle);
            end
            
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
            wPlotData = 0; %#ok<NASGU>
            
            for i=1:length(iValues)
                
                iThis.mSetHold();
                
                wPlotData = iValues{i};
                
                if (any(size(wPlotData) == 2))
                    
                    if(isa(wPlotData,'cell'))
                        %Case where we have a matrix and a character
                        %specifier
                        wMarker = wPlotData{2};
                        wPlotData = wPlotData{1};
                    else
                        wMarker = '-';
                    end
                    
                    if (size(wPlotData,1) == 2)
                        iFuncHandle(wPlotData(1,:),wPlotData(2,:),wMarker);
                    else
                        iFuncHandle(wPlotData(:,1),wPlotData(:,2),wMarker);
                    end
                    
                elseif (size(wPlotData,2) == 1)
                    iFuncHandle(wPlotData(1));
                end
                
            end
            
            if(length(varargin) > 3)
                legend(varargin{4});
                iThis.mProcessLabels(oHandle,varargin{1:3});
            else
                iThis.mProcessLabels(oHandle,varargin{:});
            end
            
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

