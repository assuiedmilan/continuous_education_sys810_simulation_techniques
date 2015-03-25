addpath(genpath('..\..'));

%Parametres
close all;
clear all; %#ok<CLSCR>
clc;

wSampleTime=0.1;
wSimulationTime=10;
wMaxStep=wSampleTime/1000;

wInputSignal=ones(1,wSimulationTime/wSampleTime);

wContinuousSystemNum = [100,0];
wContinuousSystemDen = [1,11,30,200];

wAdamsBashforthNum = [3,-1];
wAdamsBashforthDen = [2,-2,0];

%Objects initialization
wSystem = Discretizer(wSampleTime,...
    wContinuousSystemNum,...
    wContinuousSystemDen);

wPloter = Ploter([0 0 8 5],[8 5]);

%Stability study
wABStabilityHandle = wPloter.mDrawStabilityRegion('Adam-Brashforth second order',wAdamsBashforthNum,wAdamsBashforthDen);
wLambda = pole(wSystem.mGetTf);
wSampleTime=0.01:0.01:0.1;

set(0,'currentfigure',wABStabilityHandle);
for k=1:length(wSampleTime)
    
    wReal = [];
    wImag = [];
    
    for h=1:length(wLambda)
        
        wReal = [wReal,wSampleTime(k)*real(wLambda(h))];
        wImag = [wImag,wSampleTime(k)*imag(wLambda(h))];
    end
    
    hold all;
    scatter(wReal,wImag);
    legend(get(legend(gca),'String'),num2str(wSampleTime(k)));
    
end


[A,B,C,D] = wSystem.mGetStateSpaceMatrix('observable');

%Parametrage simulation
model='adamsFamillyModel';
load_system(model)
tic

wSaveFileName     = 'Y';

set_param(model,'StopFcn','save(wSaveFileName,wSaveFileName)');
set_param(strcat(model,'/Output'),'VariableName',wSaveFileName);

%Parametrage simulation continue
set_param(strcat(model,'/Continuous'),'Numerator','wContinuousSystemNum');
set_param(strcat(model,'/Continuous'),'Denominator','wContinuousSystemDen');

set_param(strcat(model,'/Observable continuous model'),'a0','A(size(A,1),1)');
set_param(strcat(model,'/Observable continuous model'),'a1','A(size(A,1),2)');
set_param(strcat(model,'/Observable continuous model'),'a2','A(size(A,1),3)');

set_param(strcat(model,'/Observable continuous model'),'b0','B(3,1)');
set_param(strcat(model,'/Observable continuous model'),'b1','B(2,1)');
set_param(strcat(model,'/Observable continuous model'),'b2','B(1,1)');

set_param(strcat(model,'/Observable continuous model'),'b3','D(1,1)');

%Parametrage AB_2
set_param(strcat(model,'/AB_2'),'a0','A(size(A,1),1)');
set_param(strcat(model,'/AB_2'),'a1','A(size(A,1),2)');
set_param(strcat(model,'/AB_2'),'a2','A(size(A,1),3)');

set_param(strcat(model,'/AB_2'),'b0','B(3,1)');
set_param(strcat(model,'/AB_2'),'b1','B(2,1)');
set_param(strcat(model,'/AB_2'),'b2','B(1,1)');

set_param(strcat(model,'/AB_2'),'b3','D(1,1)');

set_param(strcat(model,'/AB_2'),'T','wSampleTime');

set_param(strcat(model,'/AB_2'),'HzNum','wAdamsBashforthNum');
set_param(strcat(model,'/AB_2'),'HzDen','wAdamsBashforthDen');

set_param(model, 'StopTime', 'wSimulationTime');

set_param(model, 'MaxStep', 'wMaxStep');

myopts=simset('SrcWorkspace','current','DstWorkspace','current');

sim(model,wSimulationTime,myopts);
while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
end

t_sim = toc;
fprintf('\nTemps de simulation => %3.3g s\n',t_sim)

%Post traitement
load(wSaveFileName);
wStruct = eval('wSaveFileName');

%Plots

wPloter.mDrawTimeseriesPlot([Y.Continuous_signal,Y.State_space_block,Y.Observable_continuous],...
    'Open Loop Response, continuous simulation','Time (s)','Step Response');

wPloter.mDrawTimeseriesPlot([Y.Observable_continuous,Y.Observable_Adams_Branshforth],...
    'Open Loop Response, Adams Branshforth','Time (s)','Step Response','stairs');


