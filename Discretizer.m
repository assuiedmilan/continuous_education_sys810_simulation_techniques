classdef Discretizer < handle
    
    properties (SetAccess = private, GetAccess = private)
        
        mPlotMarkersList = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
        
        mT;
        mContinuousTf;
        
        %Matrices de transfer:
        mQHalijak;
        mQBoxer;
        
    end
    
    methods %Public
        
        %Constructors
        function oInstance = Discretizer(iT, varargin)
            if nargin == 2
                %Case where a TF is entered
                oInstance.mContinuousTf = varargin{1};
            elseif nargin == 3
                %Case where a NUM/DEN is entered
                oInstance.mContinuousTf = tf(varargin{1},varargin{2});
            elseif nargin == 5
                %Case where a state matrix system is entered
                [wNum,wDen] = ss2tf(varargin{1},varargin{2},varargin{3},varargin{4});
                oInstance.mContinuousTf = tf(wNum,wDen);
            else
                oInstance.mContinuousTf = 0;
                h = errordlg('Invalid constructor');
                waitfor(h);
            end
            
            oInstance.mSetSampleTime(iT);
            
        end
        
        %Public methods
        
        %Accessors
        
        function mSetSampleTime(iThis,iSampleTime)
            iThis.mT = iSampleTime;
            mUpdateMatrix(iThis);
        end
        
        function oSampleTime = mGetSampleTime(iThis)
            oSampleTime = iThis.mT;
        end
        
        function oMatrix = mGetHalijakMatrix(iThis,iRank)
            oMatrix = iThis.mQHalijak{iRank};
        end
        
        function oMatrix = mGetBoxerThalerMatrix(iThis,iRank)
            oMatrix = iThis.mQBoxer{iRank};
        end
        
        function oTf = mGetTf(iThis)
            oTf = iThis.mContinuousTf;
        end
        
        %Compute discrete TF
        function varargout = mGetDiscreteTf(iThis,iType)
            
            wDiscreteTf = mProcessTf(iThis,iThis.mContinuousTf,iType);
            
            if (nargout == 0) || (nargout == 1)
                varargout{1} = wDiscreteTf;
            elseif nargout == 2
                [wHnum,wHden] = tfdata(wDiscreteTf,'v');
                varargout{1} = wHnum./wHden(1);
                varargout{2} = wHden./wHden(1);
            else
                h = errordlg('Invalid number of output arguments, valids outputs are (Tf) or (NumTf,DenTf)');
                waitfor(h)
            end
            
        end
        
        %Get poles for the specified model of system
        function oSystemPoles = mGetPoles(iThis,iType)
            
            wSystem = mGetDiscreteTf(iThis,iType);
            
            oSystemPoles = pole(wSystem);
            
        end
        
        %Compute closed loop TF.
        function varargout = mGetClosedLoop(iThis,iFeedBackTf,iType)
            
            wDiscreteTf = mProcessTf(iThis,iThis.mContinuousTf,iType);
            wDiscreteFeedBackTf = mProcessTf(iThis,iFeedBackTf,iType);
            
            wCLTF = feedback(wDiscreteTf,wDiscreteFeedBackTf);
            
            if (nargout == 0) || (nargout == 1)
                varargout{1} = wCLTF;
            elseif nargout == 2
                [wHnum,wHden] = tfdata(wCLTF,'v');
                varargout{1} = wHnum;
                varargout{2} = wHden;
            else
                h = errordlg('Invalid number of output arguments, valids outputs are (Tf) or (NumTf,DenTf)');
                waitfor(h)
            end
            
        end
        
        %Compute discrete TF and apply retard
        function varargout = mGetRetardedDiscreteTf(iThis,iType,iRetard)
            
            wHnum    = 0; %#ok<NASGU>
            wHden    = 0; %#ok<NASGU>
            wRetard  = ones(1,iRetard);
            
            [wHnum,wHden] = mGetDiscreteTf(iThis,iType);
            wHden         = [wRetard,wHden];
            
            if (nargout == 0) || (nargout == 1)
                varargout{1} = tf(wHnum,wHden,iThis.mT);
            elseif nargout == 2
                varargout{1} = wHnum;
                varargout{2} = wHden;
            else
                h = errordlg('Invalid number of output arguments, valids outputs are (Tf) or (NumTf,DenTf)');
                waitfor(h)
            end
        end
        
        %Compute discrete TF and apply retard
        function [A,B,C,D] = mGetStateSpaceMatrix(iThis,iType)
            
            switch (iType)
                case 'observable'
                    
                    [A,B,C,D] = iThis.mProcessObservableState();
                    
                case 'commandable'
                    
                    [A,B,C,D] = iThis.mProcessCommandableState();
                    
                otherwise
                    h = errordlg('Invalid Type, returning null matrixes');
                    waitfor(h)
                    A = 0;
                    B = 0;
                    C = 0;
                    D = 0;
            end
        end
        %Execute recursion equation with the specified input, on the
        %specified typed discrete function.
        
        function Y = mComputeRecursion(iThis,U,iType)
            
            [wNum,wDen] = iThis.mGetDiscreteTf(iType);
            Y = iThis.mProcessRecursion(wNum,wDen,U);
            
        end
        
        %Stability study
        function [oStabRegionHdl,oStabCircleHdl] = mComputeStabilityRegion(iThis,iTitle,varargin)
            
            oStabRegionHdl = 0;
            oStabCircleHdl = 0;
            
            if isempty(varargin)
                h = errordlg('Error, this methods needs specific arguments, enter either a Tf for a single integrator, two Tf for a predictor-corrector, or an double for a Runge-Kutta method.');
                waitfor(h);
                return;
            end
            
            if (length(varargin) == 1 && isa(varargin{1},'tf'))
                
                %Single integrator case
                wIntegrator = Discretizer(1,varargin{:});
                [wTauCoefficients,wLCoefficients] = wIntegrator.mGetDiscreteTf('continuous');
                
                oStabRegionHdl = iThis.mProcessStabilityRegion(iTitle,wTauCoefficients,wLCoefficients);
                oStabCircleHdl = iThis.mProcessStabilityCircle(iTitle,wTauCoefficients,wLCoefficients);
                
            elseif (length(varargin) == 1 && isa(varargin{1},'double'))
                
                %Runge-Kutta
                oStabRegionHdl = iThis.mProcessRungeKuttaStabilityRegion(iTitle,varargin{1});
                
            elseif(length(varargin) == 2 && isa(varargin{1},'tf') && isa(varargin{2},'tf'))
                
                %Predictor corrector case
                wPredictor = Discretizer(1,varargin{1});
                wCorrector = Discretizer(1,varargin{2});
                
                [wTauPredictor,wRhoPredictor] = wPredictor.mGetDiscreteTf('continuous');
                [wTauCorrector,wRhoCorrector] = wCorrector.mGetDiscreteTf('continuous');
                
                oStabRegionHdl = iThis.mProcessPredictorCorrectorStabilityRegion(iTitle,wTauPredictor,wRhoPredictor,wTauCorrector,wRhoCorrector);
                oStabCircleHdl = 0;
                
            elseif(length(varargin) == 2 && isa(varargin{1},'double') && isa(varargin{2},'double'))
                %Runge-Kutta with a specified stiffed adapter
                oStabRegionHdl = iThis.mProcessRungeKuttaStabilityRegion(iTitle,varargin{1},varargin{2});
                
            end
            
        end
        
    end %Public methods
    
    %Private methods
    methods (Access = private)
        
        %Updates all matrixes when requested
        function mUpdateMatrix(iThis)
            mUpdateBoxerThalerMatrix(iThis);
            mUpdateHalijakMatrix(iThis);
        end
        
        %Halijak substitution matrix
        function mUpdateHalijakMatrix(iThis)
            
            mQ1=...
                [iThis.mT,  0;...
                1, -1];
            
            mQ2=...
                [0, iThis.mT^2, 0;...
                iThis.mT, -iThis.mT,  0;...
                1,  -2,  1];
            
            mQ3=...
                [0,  iThis.mT^3/2, iThis.mT^3/2, 0;...
                0,   iThis.mT^2,  -iThis.mT^2,   0;...
                iThis.mT, -2*iThis.mT,   iThis.mT,     0;...
                1,  -3,      3,     -1];
            
            mQ4 =...
                [0,  iThis.mT^4/4,  2*iThis.mT^4/4,  iThis.mT^4/4 , 0;...
                0,   iThis.mT^3/2,  0       , -iThis.mT^3/2 , 0;...
                0,   iThis.mT^2  , -2*iThis.mT^2  ,  iThis.mT^2   , 0;...
                iThis.mT, -3*iThis.mT  ,  3*iThis.mT    , -iThis.mT     , 0;...
                1,  -4    ,   6       , -4      , 1];
            
            iThis.mQHalijak = {mQ1, mQ2, mQ3, mQ4};
            
        end
        
        %Boxer Thalor substitution matrix
        function mUpdateBoxerThalerMatrix(iThis)
            
            mQ1=...
                [iThis.mT/2,  iThis.mT/2;...
                1         , -1];
            
            mQ2=...
                [iThis.mT^2/12, 10*iThis.mT^2/12, iThis.mT^2/12;...
                iThis.mT/2   , 0              , -iThis.mT/2;...
                1,            -2              ,  1];
            
            mQ3=...
                [0           ,    iThis.mT^3/2 ,    iThis.mT^3/2 ,  0;...
                iThis.mT^2/12,  9*iThis.mT^2/12, -9*iThis.mT^2/12, -iThis.mT^2/12;...
                iThis.mT/2   , -1*iThis.mT/2   , -1*iThis.mT/2   ,  iThis.mT/2;...
                1            , -3             ,   3              , -1];
            
            mQ4 =...
                [-iThis.mT^4/720,  124*iThis.mT^4/720,  474*iThis.mT^4/720,  124*iThis.mT^4/720, -iThis.mT^4/720;...
                0               ,   iThis.mT^3/2     ,  0                 , -iThis.mT^3/2      , 0;...
                iThis.mT^2/12   ,   8*iThis.mT^2/12  , -18*iThis.mT^2/12  ,  8*iThis.mT^2/12   , iThis.mT^2/12;...
                iThis.mT/2      , -iThis.mT          , 0                  ,  iThis.mT          , -iThis.mT/2;...
                1               ,  -4                ,   6                , -4                 , 1];
            
            iThis.mQBoxer = {mQ1, mQ2, mQ3, mQ4};
            
        end
        
        %Execute recursion equation with the specified input.
        function Y = mProcessRecursion(iThis,num,den,U)
            
            iThis; %#ok<VUNUS>
            
            wTrimedDen = den(find(den,1):size(den,2));
            wTrimedNum = num(find(num,1):size(num,2));
            
            %How many iterations of Y are not computable.
            %If the system is implicit, then wUdelta = 0
            wUdelta = size(wTrimedDen,2)-size(wTrimedNum,2);
            if (wUdelta < 0)
                h = errordlg('Error, non-causal system');
                waitfor(h);
                return;
            end
            
            Y = zeros(1,wUdelta);
            
            for i = size(Y,2)+1:size(U,2)
                
                wY = 0;
                wU = 0;
                
                %Building input sum according to matlab 1-based indexing
                for k = 1:size(wTrimedNum,2)
                    if(i-k+1-wUdelta > 0)
                        wU = wU + wTrimedNum(k)*U(i-k+1-wUdelta);
                    else
                        wU = wU + 0;
                    end
                end
                
                %Building output sum according to matlab 1-based indexing
                for k = 1:size(wTrimedDen,2)-1
                    if((i-k > 0) && (i-k <= size(Y,2)))
                        wY = wY + wTrimedDen(k+1)*Y(i-k);
                    else
                        wY = wY + 0;
                    end
                end
                
                %Sum must be pondered by the highest numerator coefficient
                Y(i) = 1/wTrimedDen(1) * (-wY + wU);
                
            end
        end
        
        %Process TF conversions.
        function oTf = mProcessTf(iThis,iTf,iType)
            
            wTf      = iTf;
            wHnum    = 0; %#ok<NASGU>
            wHden    = 0; %#ok<NASGU>
            
            if (strcmp(get(iTf,'Variable'),'s'))
                switch (iType)
                    case 'zoh'
                        
                        wTf = c2d(iTf,iThis.mT,'zoh');
                        
                    case 'tutsin'
                        
                        wTf = c2d(iTf,iThis.mT,'tutsin');
                        
                    case {'halijak','boxerThaler'}
                        
                        [wHnum,wHden] = tfdata(iTf,'v');
                        
                        if(strcmp(iType,'halijak'))
                            wH = [fliplr(wHnum);fliplr(wHden)]*iThis.mGetHalijakMatrix(size(wHden,2)-1);
                        else
                            wH = [fliplr(wHnum);fliplr(wHden)]*iThis.mGetBoxerThalerMatrix(size(wHden,2)-1);
                        end
                        
                        wTf = tf(wH(1,:),wH(2,:),iThis.mT);
                        
                    case 'continuous'
                        %Return input
                        
                    otherwise
                        warning('mProcessTf returning input (continuous) TF')
                end
                
            else
                switch (iType)
                    case 'zoh'
                        
                        wTf = d2d(iTf,iThis.mT,'zoh');
                        
                    case 'tutsin'
                        
                        wTf = d2d(iTf,iThis.mT,'tutsin');
                        
                    otherwise
                        
                        warning('mProcessTf returning input (discrete) TF')
                end
                
            end
            
            oTf = wTf;
            
        end
        
        %Process state space
        function [A,B,C,D] = mProcessObservableState(iThis)
            
            wOneBaseBias = 1;
            
            [wNum,wDen,wOrder] = iThis.mGetMatlabOrderedCoefficients();
            
            A = iThis.mProcessCannonicalAMatrix();
            
            %Compute beta values
            wBetaMatrix = zeros(wOrder+1,1);
            for k=wOrder:-1:0
                s= 0;
                for i=1:wOrder-k
                    s = s + wDen(wOrder-i+wOneBaseBias)*wBetaMatrix(k+i+wOneBaseBias);
                end
                
                wBetaMatrix(k+wOneBaseBias) = wNum(k+wOneBaseBias) - s;
            end
            
            %Extract B & D matrixes from beta values. Note that beta matrix
            %is ordered in matlab order (aka: [b(0);b(1);...b(n-1);b(n)])
            %B needs to be inverted to be in the following order:
            %[b(n-1),b(n-2),...b(0)]
            B = wBetaMatrix(length(wBetaMatrix)-1:-1:1);
            D = wBetaMatrix(length(wBetaMatrix));
            
            %Compute C matrix.
            C = [1,zeros(1,wOrder-1)];
        end
        
        function [A,B,C,D] = mProcessCommandableState(iThis)
            
            wOneBaseBias = 1;
            
            [wNum,wDen,wOrder] = iThis.mGetMatlabOrderedCoefficients(); %#ok<ASGLU>
            
            A = iThis.mProcessCannonicalAMatrix();
            
            %Compute C values
            C = zeros(1,wOrder);
            for i=0:wOrder-1
                C(1,i+wOneBaseBias) = wNum(i+wOneBaseBias);
            end
            
            %Compute B and D
            B = [zeros(1,wOrder-1),1]';
            D = 0;
            
        end
        
        function [A] = mProcessCannonicalAMatrix(iThis)
            
            wOneBaseBias = 1;
            
            [wNum,wDen,wOrder] = iThis.mGetMatlabOrderedCoefficients(); %#ok<ASGLU>
            
            %Compute A matrix
            A = diag(ones(1,wOrder-1),1);
            for i=0:wOrder-1
                A(wOrder,i+wOneBaseBias) = -wDen(i+wOneBaseBias);
            end
            
        end
        
        function [oNum,oDen,oOrder] = mGetMatlabOrderedCoefficients(iThis)
            
            %Get transfert function polynomials and set them in a proper
            %matlab order, opposed to polynomial order (aka [a(0),a(1),...,a(n-1),a(n)])
            
            %Also normalize the coefficients
            [oNum,oDen] = tfdata(iThis.mGetTf());
            
            oDen = fliplr(oDen{1});
            oNum = fliplr(oNum{1});
            
            oDen = oDen/oDen(length(oDen));
            oNum = oNum/oDen(length(oDen));
            
            oOrder = length(oDen)-1;
            
        end
        
        function mAddPolesToStabilityRegion(iThis,iFigureHandle,iPloter)
            
            wSystemPoles = iThis.mGetPoles('continuous');
            wSampleTime = iThis.mGetSampleTime();
            
            for k=1:length(wSystemPoles)
                
                hold all;
                wPlotMarkersListIndex = mod(k,length(iThis.mPlotMarkersList));
                plot(real(wSystemPoles(k))*wSampleTime,imag(wSystemPoles(k))*wSampleTime,iThis.mPlotMarkersList(wPlotMarkersListIndex));
                legend(get(legend(gca),'String'),['Pole ',num2str(wSystemPoles(k))]);
                
            end
            
            iPloter.mProcessSaveDraw(iFigureHandle);
        end
        
        
        function oHandle = mProcessRungeKuttaStabilityRegion(iThis,iTitle,iOrder,varargin)
            
            wMatlabIndexBias = 1;
            
                       
            wZvalues = (exp(1i*(0:0.01:2*pi)))';
            
            wPloter = Ploter([0 0 5 5],[5 5]);
            wPloter.mSetSaveAfterDraw(false);
            
            if(not(isempty(varargin)))
                wStiffAdapterMatrix = varargin{1};
                wNumberOfPlots = length(wStiffAdapterMatrix);
                wPonderingCoefficients = [1];
            else
                if(exist('wStiffAdapterMatrix')) %#ok<EXIST>
                    clear wStiffAdapterMatrix
                end
                wNumberOfPlots = 1;
                wPonderingCoefficients = [1];
            end
            
            wPlotDataCell = cell(1,wNumberOfPlots);
            for wCurrentPlot = 1:wNumberOfPlots
                
                clear wStabilityPolynom wStabilityValues wPlotData
                wStabilityPolynom = zeros(1,iOrder+1);
                wStabilityValues = zeros(iOrder*length(wZvalues),1);
                wPlotData = zeros(length(wZvalues)*iOrder*length(wPonderingCoefficients),2);
                            
                for k=1:iOrder
                    if (exist('wStiffAdapterMatrix') && k==iOrder)                         %#ok<EXIST>
                        wStabilityPolynom(k+wMatlabIndexBias) = 1/wStiffAdapterMatrix(wCurrentPlot);
                    else
                        wStabilityPolynom(k+wMatlabIndexBias) = 1/factorial(k);
                    end
                end
                
                %Correctly order coefficients for matlab
                wStabilityPolynom = fliplr(wStabilityPolynom);
                
                for l = 1:length(wPonderingCoefficients)
                    for k = 1:length(wZvalues)
                        
                        wStabilityPolynom(iOrder+1) = 1-wPonderingCoefficients(l)*wZvalues(k);
                        wRoots = roots(wStabilityPolynom);
                        
                        for h=0:iOrder-1
                            wStabilityValues(iOrder*k-h,l) = wRoots(h+wMatlabIndexBias);
                        end
                        
                    end
                    
                    wPlotData((l-1)*length(wZvalues)*iOrder + wMatlabIndexBias:l*length(wZvalues)*iOrder,1:2) = [real(wStabilityValues(:,l)),imag(wStabilityValues(:,l))];
                end
                
                wPlotDataCell{1,wCurrentPlot} = {wPlotData,'.'};
            end
            
            oHandle = wPloter.mDrawStandardPlot(wPlotDataCell...
                ,'plot'...
                ,['Stability Region ',iTitle]...
                ,'Real axis'...
                ,'Imaginary axis'...
                ,'Stability region');
            
            iThis.mAddPolesToStabilityRegion(oHandle,wPloter);
        end
        
        function oHandle = mProcessPredictorCorrectorStabilityRegion(iThis,iTitle,iTauPredictor,iRhoPredictor,iTauCorrector,iRhoCorrector)
            
            wZvalues = [(exp(1i*(0:.01:2*pi)))',0.9*(exp(1i*(0:.01:2*pi)))'];
            wPloter = Ploter([0 0 5 5],[5 5]);
            wPloter.mSetSaveAfterDraw(false);
            
            wBetaK = iTauCorrector(1);
            
            wStabilityRoots = zeros(2*size(wZvalues,1),size(wZvalues,2));
            for l=1:size(wZvalues,2)
                for k=1:size(wZvalues,1)
                    
                    wP1 = -wBetaK * polyval(iTauPredictor,wZvalues(k,l));
                    wP2 = wBetaK * polyval(iRhoPredictor,wZvalues(k,l)) - polyval(iTauCorrector,wZvalues(k,l));
                    wP3 = polyval(iRhoCorrector,wZvalues(k,l));
                    
                    wRoots = roots([wP1,wP2,wP3]);
                    wStabilityRoots(2*k-1,l) = wRoots(1);
                    wStabilityRoots(2*k,l) = wRoots(2);
                    
                end
            end
            
            oHandle= wPloter.mDrawStandardPlot({{[real(wStabilityRoots(:,1)),imag(wStabilityRoots(:,1))],'.'}...
                ,{[real(wStabilityRoots(:,2)),imag(wStabilityRoots(:,2))],'.'}}...
                ,'plot'...
                ,['Stability Region ',iTitle,strrep(num2str(iThis.mGetSampleTime()*1000),'.',''),'ms']...
                ,'Real axis'...
                ,'Imaginary axis'...
                ,{'z=e(i.wt)','z=0.9*e(i.wt)'});
            
            iThis.mAddPolesToStabilityRegion(oHandle,wPloter);
            
        end
        
        function oHandle = mProcessStabilityRegion(iThis,iTitle,iTauCoefficients,iLCoefficients)
            
            wZvalues = (exp(1i*(0:.01:2*pi)))';
            wPloter = Ploter([0 0 5 5],[5 5]);
            wPloter.mSetSaveAfterDraw(false);
            
            wStabilityValues = zeros(1,length(wZvalues));
            for k = 1:length(wZvalues)
                
                wStabilityValues(k) = polyval(iLCoefficients,wZvalues(k)) ./ polyval(iTauCoefficients,wZvalues(k));
                
            end
            
            oHandle = wPloter.mDrawStandardPlot({[real(wStabilityValues);imag(wStabilityValues)]}...
                ,'plot'...
                ,['Stability Region ',iTitle,strrep(num2str(iThis.mGetSampleTime()*1000),'.',''),'ms']...
                ,'Real axis'...
                ,'Imaginary axis'...
                ,'Stability region');
            
            iThis.mAddPolesToStabilityRegion(oHandle,wPloter);
        end
        
        function oHandle = mProcessStabilityCircle(iThis,iTitle,iTauCoefficients,iLCoefficients)
            
            %Unit circle param
            Wteta=0:0.001:2*pi;
            wX=1*cos(Wteta);
            wY=1*sin(Wteta);
            wSampleTime = iThis.mGetSampleTime();
            
            wLegendString = {'Stability Circle'};
            
            wPloter = Ploter([0 0 5 5],[5 5]);
            wSystemPoles = iThis.mGetPoles('continuous');
            
            wExactSystemPoles = exp(wSystemPoles*wSampleTime);
            
            wClosedLoopPoles = cell(1,2*length(wExactSystemPoles));
            
            for k=1:1:length(wExactSystemPoles)
                
                %Storing discret poles
                wClosedLoopRoots = roots(iLCoefficients-wExactSystemPoles(k)*wSampleTime*iTauCoefficients);
                
                wPlotMarkersListIndex = mod(k,length(iThis.mPlotMarkersList));
                
                wRealPart = zeros(1,length(wClosedLoopRoots));
                wImagPart = zeros(1,length(wClosedLoopRoots));
                
                for j=1:length(wClosedLoopRoots)
                    wRealPart(j) = real(wClosedLoopRoots(j));
                    wImagPart(j) = imag(wClosedLoopRoots(j));
                end
                
                wClosedLoopPoles{2*k-1} = {[wRealPart;wImagPart],iThis.mPlotMarkersList(wPlotMarkersListIndex)};
                wLegendString{2*k} = ['Pole: ',num2str(wSystemPoles(k))];
                
                %Storing exact pole on next index
                
                wPlotMarkersListIndex = mod(k+1,length(iThis.mPlotMarkersList));
                wClosedLoopPoles{2*k} = {[real(wExactSystemPoles(k));imag(wExactSystemPoles(k))],iThis.mPlotMarkersList(wPlotMarkersListIndex)};
                wLegendString{2*k+1} = ['Exact Pole: ',num2str(wSystemPoles(k))];
                
            end
            
            oHandle = wPloter.mDrawStandardPlot([[wX;wY],wClosedLoopPoles]...
                ,'plot'...
                ,['Stability Circle ',iTitle,strrep(num2str(wSampleTime*1000),'.',''),'ms']...
                ,'Real axis'...
                ,'Imaginary axis'...
                ,wLegendString);
        end
        
    end %Private methods
    
end %Class

