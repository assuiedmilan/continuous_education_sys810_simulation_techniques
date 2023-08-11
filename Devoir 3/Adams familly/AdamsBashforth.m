addpath(genpath('..\..'));
dbstop if error

% ************************************************************ %
% ****************      INITIALIZATION      ****************** %
% ************************************************************ %

close all;
clear all; %#ok<CLSCR>
clc;

wSampleTimes=[0.2,0.1,0.09,0.05];
wSimulationTime=10;

wContinuousSystemNum = [100,0];
wContinuousSystemDen = [1,11,30,200];

wAdamsBashforthNum = [3,-1];
wAdamsBashforthDen = [2,-2,0];

wPloter = Ploter([0 0 5 5],[5 5]);

for i=1:length(wSampleTimes)
    
    close all;
    java.lang.Runtime.getRuntime.freeMemory;
    
    clearvars -except wSampleTimes wSimulationTime wContinuousSystemNum...
    wContinuousSystemDen wAdamsBashforthNum wAdamsBashforthDen wPloter i
    
    wSampleTime = wSampleTimes(i);
    wMaxStep=wSampleTime/1000;

    %Objects initialization
    wSystem = Discretizer(wSampleTime,...
        wContinuousSystemNum,...
        wContinuousSystemDen);
    
    [Ao,Bo,Co,Do] = wSystem.mGetStateSpaceMatrix('observable');
    [Ac,Bc,Cc,Dc] = wSystem.mGetStateSpaceMatrix('commandable');
    
    %Stability study
    wSystem.mComputeStabilityRegion('Adam-Brashforth second order ',tf(wAdamsBashforthNum,wAdamsBashforthDen));
    
    % ************************************************************ %
    % **********************  SIMULATION ************************* %
    % ************************************************************ %
    adamsFamillyModelSetup;
	
    myopts=simset('SrcWorkspace','current','DstWorkspace','current');
    
    sim(model,wSimulationTime,myopts);
    while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
    end
    
    t_sim = toc;
    fprintf('\nTemps de simulation => %3.3g s\n',t_sim)
    
    %Drawing and saving Plots
    
    wPloter.mDrawTimeseriesPlot([SimOutput.Continuous_signal...
        ,SimOutput.Commandable_continuous...
        ,SimOutput.Observable_continuous]...
        ,['Continuous simulation ',strrep(num2str(wSampleTime*1000),'.',''),'ms']...
        ,'Time (s)'...
        ,'Step Response');
    
    wPloter.mDrawTimeseriesPlot([SimOutput.Observable_continuous...
        ,SimOutput.Observable_Adams_Branshforth]...
        ,['Adams Branshforth ',strrep(num2str(wSampleTime*1000),'.',''),'ms']...
        ,'Time (s)'...
        ,'Step Response'...
        ,'stairs');
    
    close all;    
end

%Saving Model
wPloter.mProcessSaveModel(model);
wPloter.mProcessSaveModel(wObsBlock);
wPloter.mProcessSaveModel(wComBlock);
wPloter.mProcessSaveModel(wAB2Block);




