addpath(genpath('..\..'));

% ************************************************************ %
% ****************      INITIALIZATION      ****************** %
% ************************************************************ %

close all;
clear all; %#ok<CLSCR>
clc;

wSampleTimes=[0.5,0.10,0.09,0.01];
wSimulationTime=10;

wContinuousSystemNum = [100,0];
wContinuousSystemDen = [1,11,30,200];

wAdamsBashforthNum = [3,-1];
wAdamsBashforthDen = [2,-2,0];

wPloter = Ploter([0 0 8 5],[8 5]);

for i=1:length(wSampleTimes)
    
    wSampleTime = wSampleTimes(i);
    wMaxStep=wSampleTime/1000;

    %Objects initialization
    wSystem = Discretizer(wSampleTime,...
        wContinuousSystemNum,...
        wContinuousSystemDen);
    [Ao,Bo,Co,Do] = wSystem.mGetStateSpaceMatrix('observable');
    [Ac,Bc,Cc,Dc] = wSystem.mGetStateSpaceMatrix('commandable');
    
    %Stability study
    wSystem.mComputeStabilityRegion('Adam-Brashforth second order',wSampleTime,wAdamsBashforthNum,wAdamsBashforthDen);
    
    % ************************************************************ %
    % ****************  SIMULINK INITIALIZATION ****************** %
    % ************************************************************ %
    
    model='adamsFamillyModel';
    load_system(model)
    tic
    
    wSaveFileName     = 'SimOutput';
    wContBlock = strcat(model,'/Continuous');
    wObsBlock = strcat(model,'/Observable continuous');
    wComBlock = strcat(model,'/Commandable continuous');
    wAB2Block = strcat(model,'/AB_2');
    
    set_param(model,'StopFcn','save(wSaveFileName,wSaveFileName)');
    set_param(strcat(model,'/Output'),'VariableName',wSaveFileName);
    
    %Parametrage simulation continue
    set_param(wContBlock,'Numerator','wContinuousSystemNum');
    set_param(wContBlock,'Denominator','wContinuousSystemDen');
    
    set_param(wObsBlock,'a0','Ao(size(Ao,1),1)');
    set_param(wObsBlock,'a1','Ao(size(Ao,1),2)');
    set_param(wObsBlock,'a2','Ao(size(Ao,1),3)');
    
    set_param(wObsBlock,'b0','Bo(3,1)');
    set_param(wObsBlock,'b1','Bo(2,1)');
    set_param(wObsBlock,'b2','Bo(1,1)');
    
    set_param(wObsBlock,'b3','Do(1,1)');
    
    set_param(wComBlock,'a0','Ac(size(Ac,1),1)');
    set_param(wComBlock,'a1','Ac(size(Ac,1),2)');
    set_param(wComBlock,'a2','Ac(size(Ac,1),3)');
    
    set_param(wComBlock,'b0','Cc(1,1)');
    set_param(wComBlock,'b1','Cc(1,2)');
    set_param(wComBlock,'b2','Cc(1,3)');
    
    %Parametrage AB_2
    set_param(wAB2Block,'a0','Ao(size(Ao,1),1)');
    set_param(wAB2Block,'a1','Ao(size(Ao,1),2)');
    set_param(wAB2Block,'a2','Ao(size(Ao,1),3)');
    
    set_param(wAB2Block,'b0','Bo(3,1)');
    set_param(wAB2Block,'b1','Bo(2,1)');
    set_param(wAB2Block,'b2','Bo(1,1)');
    
    set_param(wAB2Block,'b3','Do(1,1)');
    
    set_param(wAB2Block,'T','wSampleTime');
    
    set_param(wAB2Block,'HzNum','wAdamsBashforthNum');
    set_param(wAB2Block,'HzDen','wAdamsBashforthDen');
    
    set_param(model, 'StopTime', 'wSimulationTime');
    
    set_param(model, 'MaxStep', 'wMaxStep');
    
    % ************************************************************ %
    % **********************  SIMULATION ************************* %
    % ************************************************************ %
    
    myopts=simset('SrcWorkspace','current','DstWorkspace','current');
    
    sim(model,wSimulationTime,myopts);
    while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
    end
    
    t_sim = toc;
    fprintf('\nTemps de simulation => %3.3g s\n',t_sim)
    
    % ************************************************************ %
    % ****************  PREDICTION - CORRECTION ****************** %
    % ************************************************************ %
    
    %Parameters
    t=linspace(0,wSimulationTime,wSimulationTime*(1/wSampleTime));
    
    %Initial conditions
    f1c(1) = 0;   f2c(1) = 0; f3c(1) = 1;
    x1c(1) = 0;   x2c(1) = 0; x3c(1) = 0;
    
    for n = 0:wSimulationTime/wSampleTime-2
        
        if (n == 0)
            x1p(n+2) = x1c(n+1) + (wSampleTime/2)*(f1c(n+1)); %#ok<*SAGROW>
            x2p(n+2) = x2c(n+1) + (wSampleTime/2)*(f2c(n+1));
            x3p(n+2) = x3c(n+1) + (wSampleTime/2)*(f3c(n+1));
        else
            x1p(n+2) = x1c(n+1) + (wSampleTime/2)*(f1c(n+1)-f1c(n));
            x2p(n+2) = x2c(n+1) + (wSampleTime/2)*(f2c(n+1)-f2c(n));
            x3p(n+2) = x3c(n+1) + (wSampleTime/2)*(f3c(n+1)-f3c(n));
        end
        
        f1p(n+2) = x2p(n+2);
        f2p(n+2) = x3p(n+2);
        f3p(n+2) = -200*x1p(n+2) - 30*x2p(n+2) -11*x3p(n+2) + 1;
        
        x1c(n+2) = x1c(n+1) + (wSampleTime/2)*(f1p(n+2)+ f1c(n+1));
        x2c(n+2) = x2c(n+1) + (wSampleTime/2)*(f2p(n+2)+ f2c(n+1));
        x3c(n+2) = x3c(n+1) + (wSampleTime/2)*(f3p(n+2)+ f3c(n+1));
        
        f1c(n+2) = x2c(n+2);
        f2c(n+2) = x3c(n+2);
        f3c(n+2) = -200*x1c(n+2) - 30*x2c(n+2) -11*x3c(n+2) + 1;
        
    end
    
    % ************************************************************ %
    % *******************  EXTRACTING FIGURES ******************** %
    % ************************************************************ %
    
    %Drawing and saving Plots
    
    wPloter.mDrawTimeseriesPlot([SimOutput.Continuous_signal...
        ,SimOutput.Commandable_continuous...
        ,SimOutput.Observable_continuous]...
        ,['Open Loop Response, continuous simulation T',strrep(num2str(wSampleTime),'.','-'),'s']...
        ,'Time (s)'...
        ,'Step Response');
    
    wPloter.mDrawTimeseriesPlot([SimOutput.Observable_continuous...
        ,SimOutput.Observable_Adams_Branshforth]...
        ,['Open Loop Response, Adams Branshforth T',strrep(num2str(wSampleTime),'.','-'),'s']...
        ,'Time (s)'...
        ,'Step Response'...
        ,'stairs');
    
    wPloter.mDrawStandardPlot({[SimOutput.Observable_continuous.Time...
        ,SimOutput.Observable_continuous.Data]...
        ,[t;100*x2c]}...
        ,'stairs'...
        ,['Open Loop Response, Prediction-Correction vs Continuous simulation T',strrep(num2str(wSampleTime),'.','-'),'s']...
        ,'Time (s)'...
        ,'Step Response'...
        ,{'Commandable continuous';'Prediction-Correction'});
    
end

%Saving Model
wPloter.mProcessSaveModel(model);
wPloter.mProcessSaveModel(wObsBlock);
wPloter.mProcessSaveModel(wComBlock);
wPloter.mProcessSaveModel(wAB2Block);




