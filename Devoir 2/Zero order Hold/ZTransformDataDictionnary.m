addpath(genpath('..'));

%Parametres
close all;
clear all; %#ok<CLSCR>
clc;

wSampleTime=0.1;
wSimulationTime=10;
wMaxStep=wSampleTime/1000;

wInputSignal=ones(1,wSimulationTime/wSampleTime);

wContinuousSystemNum = 100;
wContinuousSystemDen = [1,11,30,200];

wDiscretizer = Discretizer(wSampleTime,...
    wContinuousSystemNum,...
    wContinuousSystemDen);

% Getting recursive equation and Tf for Zero Order Hold:
wYzod = wDiscretizer.mComputeRecursion(wInputSignal,'zoh');
[wYzodN,wYzodD] = wDiscretizer.mGetDiscreteTf('zoh');

% Getting recursive equation and Tf for Tutsin:
wYtut = wDiscretizer.mComputeRecursion(wInputSignal,'tutsin');
[wYtutN,wYtutD] = wDiscretizer.mGetDiscreteTf('tutsin');

% Getting recursive equation and Tf for Halijak:
wYhal = wDiscretizer.mComputeRecursion(wInputSignal,'halijak');
[wYhalN,wYhalD] = wDiscretizer.mGetDiscreteTf('halijak');

% Getting recursive equation and Tf for Boxer Thaler:
wYbox = wDiscretizer.mComputeRecursion(wInputSignal,'boxerThaler');
[wYboxN,wYboxD] = wDiscretizer.mGetDiscreteTf('boxerThaler');

%Parametrage simulation
model='zodModel';
load_system(model)
tic

wSaveFileName     = 'Ymodel';

set_param(model,'StopFcn','save(wSaveFileName,wSaveFileName)');
set_param(strcat(model,'/Output'),'VariableName',wSaveFileName);

set_param(strcat(model,'/ZOD'),'Numerator','wYzodN');
set_param(strcat(model,'/ZOD'),'Denominator','wYzodD');
set_param(strcat(model,'/ZOD'),'SampleTime','wSampleTime');

set_param(strcat(model,'/Tutsin'),'Numerator','wYtutN');
set_param(strcat(model,'/Tutsin'),'Denominator','wYtutD');
set_param(strcat(model,'/Tutsin'),'SampleTime','wSampleTime');

set_param(strcat(model,'/Halikaj'),'Numerator','wYhalN');
set_param(strcat(model,'/Halikaj'),'Denominator','wYhalD');
set_param(strcat(model,'/Halikaj'),'SampleTime','wSampleTime');

set_param(strcat(model,'/Boxer Thalor'),'Numerator','wYboxN');
set_param(strcat(model,'/Boxer Thalor'),'Denominator','wYboxD');
set_param(strcat(model,'/Boxer Thalor'),'SampleTime','wSampleTime');

set_param(strcat(model,'/Continuous'),...
    'Numerator','wContinuousSystemNum');
set_param(strcat(model,'/Continuous'),...
    'Denominator','wContinuousSystemDen');

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
h_zod=figure();
hold all
plot(Ymodel.Yzod.Time,Ymodel.Yzod.Data);
stairs(Ymodel.Yzod.Time,Ymodel.Yzod.Data);
stairs(0:wSampleTime:(size(wInputSignal,2)-1)*wSampleTime,wYzod);
legend('Simulink ZOH (plot)',...
    'Simulink ZOH (stairs)',...
    'Recursive equation ZOH');
grid minor;
xlabel('Time (s)');
ylabel('Step Response');
title('Open Loop Response to ZOH discretized function transfer');

h_tut=figure();
hold all
plot(Ymodel.Ytut.Time,Ymodel.Ytut.Data);
stairs(Ymodel.Ytut.Time,Ymodel.Ytut.Data);
stairs(0:wSampleTime:(size(wInputSignal,2)-1)*wSampleTime,wYtut);
plot(Ymodel.Ys.Time,Ymodel.Ys.Data);
legend('Simulink Tutsin (plot)',...
    'Simulink Tutsin (stairs)',...
    'Recursive equation Tutsin',...
    'Continuous system');
grid minor;
xlabel('Time (s)');
ylabel('Step Response');
title('Open Loop Response to Tutsin discretized function transfer');

h_hal=figure();
hold all
plot(Ymodel.Yhal.Time,Ymodel.Yhal.Data);
stairs(Ymodel.Yhal.Time,Ymodel.Yhal.Data);
stairs(0:wSampleTime:(size(wInputSignal,2)-1)*wSampleTime,wYhal);
plot(Ymodel.Ys.Time,Ymodel.Ys.Data);
legend('Simulink Halijak (plot)',...
    'Simulink Halijak (stairs)',...
    'Recursive equation Halijak','Continuous system');
grid minor;
xlabel('Time (s)');
ylabel('Step Response');
title('Open Loop Response to Halijak discretized function transfer');

h_box=figure();
hold all
plot(Ymodel.Ybox.Time,Ymodel.Ybox.Data);
stairs(Ymodel.Ybox.Time,Ymodel.Ybox.Data);
stairs(0:wSampleTime:(size(wInputSignal,2)-1)*wSampleTime,wYbox);
plot(Ymodel.Ys.Time,Ymodel.Ys.Data);
legend('Simulink Boxer-Thaler (plot)',...
    'Simulink Boxer-Thaler (stairs)',...
    'Recursive equation Boxer-Thaler','Continuous system');
grid minor;
xlabel('Time (s)');
ylabel('Step Response');
title('Open Loop Response to Boxer-Thaler discretized function transfer');

%Poles - zeros
h_spz = figure();
pzmap(wDiscretizer.mGetTf());
legend('Continuous system');
sgrid;
grid minor;

h_zpz  = figure();
pzmap(wDiscretizer.mGetDiscreteTf('zoh'),...
    wDiscretizer.mGetDiscreteTf('tutsin'),...
    wDiscretizer.mGetDiscreteTf('halijak'),...
    wDiscretizer.mGetDiscreteTf('boxerThaler'));
legend('Zero order hold','Tutsin','Halijak','Boxer-Thaler');
zgrid;
grid minor;

h_snyq = figure();
nyquist(wDiscretizer.mGetTf());
legend('Continuous system');
grid minor;

%Saving figures
wPaperPos = [0 0 8 5];
wPaperSize = [8 5];

set(h_zod, 'PaperPosition', wPaperPos); 
set(h_zod, 'PaperSize', wPaperSize);
saveas(h_zod, 'ZTransform ZoH', 'pdf')

set(h_tut, 'PaperPosition', wPaperPos); 
set(h_tut, 'PaperSize', wPaperSize); 
saveas(h_tut, 'ZTransform Tutsin', 'pdf') 

set(h_hal, 'PaperPosition', wPaperPos); 
set(h_hal, 'PaperSize', wPaperSize); 
saveas(h_hal, 'ZTransform Halijak', 'pdf') 

set(h_box, 'PaperPosition', wPaperPos); 
set(h_box, 'PaperSize', wPaperSize); 
saveas(h_box, 'ZTransform Boxer-Thalor', 'pdf') 

set(h_spz, 'PaperPosition', wPaperPos); 
set(h_spz, 'PaperSize', wPaperSize); 
saveas(h_spz, 'ZTransform Poles-Zeros', 'pdf') 

set(h_zpz, 'PaperPosition', wPaperPos); 
set(h_zpz, 'PaperSize', wPaperSize); 
saveas(h_zpz, 'ZTransform Poles-Zeros discretes', 'pdf') 

set(h_snyq, 'PaperPosition', wPaperPos); 
set(h_snyq, 'PaperSize', wPaperSize); 
saveas(h_snyq, 'ZTransform Nyquist', 'pdf') 

%Saving model
saveas(get_param(model,'Handle'),'ZTransform-Model.pdf');
close_system(model,false)

%Print the transfer functions on command windows for saving:
clc;
disp('Continuous transfer function')
wDiscretizer.mGetTf

disp('Tutsin discretization')
wDiscretizer.mGetDiscreteTf('tutsin')

disp('Boxer-Thalor discretization')
wDiscretizer.mGetDiscreteTf('boxerThaler')

disp('Halijak discretization')
wDiscretizer.mGetDiscreteTf('halijak')