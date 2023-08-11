addpath(genpath('..'));

%Parametres
close all;
clear all; %#ok<CLSCR>
clc;

%Parametres
wSampleTime = 0.04;
wSimulationTime=20;
wMaxStep=wSampleTime/2;

%Fonctions de transfert
wSysNum  = [0,0,9.8,-4.9,61.4];
wSysDen  = [1,0.44,-0.007,0.11,0];
wCompNum = [0,1.93,1.72,0.43,0.11];
wCompDen = [1,9.76,40.9,76.6,136];

wSystem      = Discretizer(wSampleTime,wSysNum,wSysDen);
wCompensator = Discretizer(wSampleTime,wCompNum,wCompDen);

%Tutsin:
[wSys_Tutn,wSys_Tutd]   = wSystem.mGetDiscreteTf('tutsin');
[wComp_Tutn,wComp_Tutd] = wCompensator.mGetDiscreteTf('tutsin');

%Halijak
[wSys_Haln,wSys_Hald]   = wSystem.mGetDiscreteTf('halijak');
[wComp_Haln,wComp_Hald] = wCompensator.mGetDiscreteTf('halijak');


%Parametrage simulation
model='HelicoModel';
load_system(model)
tic

wSaveFileName     = 'Ymodel';

set_param(model,'StopFcn','save(wSaveFileName,wSaveFileName)');
set_param(strcat(model,'/Output'),'VariableName',wSaveFileName);

set_param(strcat(model,'/Tutsin implicit'),'S_num','wSys_Tutn');
set_param(strcat(model,'/Tutsin implicit'),'S_den','wSys_Tutd');
set_param(strcat(model,'/Tutsin implicit'),'C_num','wComp_Tutn');
set_param(strcat(model,'/Tutsin implicit'),'C_den','wComp_Tutd');
set_param(strcat(model,'/Tutsin implicit'),'T','wSampleTime');

set_param(strcat(model,'/Tutsin explicit'),'S_num','wSys_Tutn');
set_param(strcat(model,'/Tutsin explicit'),'S_den','wSys_Tutd');
set_param(strcat(model,'/Tutsin explicit'),'C_num','wComp_Tutn');
set_param(strcat(model,'/Tutsin explicit'),'C_den','wComp_Tutd');
set_param(strcat(model,'/Tutsin explicit'),'T','wSampleTime');

set_param(strcat(model,'/Halijak'),'S_num','wSys_Haln');
set_param(strcat(model,'/Halijak'),'S_den','wSys_Hald');
set_param(strcat(model,'/Halijak'),'C_num','wComp_Tutn');
set_param(strcat(model,'/Halijak'),'C_den','wComp_Tutd');
set_param(strcat(model,'/Halijak'),'T','wSampleTime');

set_param(strcat(model,'/Continuous'),'S_num','wSysNum');
set_param(strcat(model,'/Continuous'),'S_den','wSysDen');
set_param(strcat(model,'/Continuous'),'C_num','wCompNum');
set_param(strcat(model,'/Continuous'),'C_den','wCompDen');

set_param(model, 'StopTime', 'wSimulationTime');
set_param(model, 'MaxStep', 'wMaxStep');
set_param((strcat(model,'/Output')), 'SampleTime','wSampleTime');

myopts=simset('SrcWorkspace','current','DstWorkspace','current');

sim(model,wSimulationTime,myopts);
while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
end

t_sim = toc;
fprintf('\nTemps de simulation => %3.3g s\n',t_sim)

%Post traitement
load(wSaveFileName);
wStruct = eval(wSaveFileName);

%Plots
h_zod=figure();
hold all
plot(Ymodel.Tutsin_explicit_error.Time,...
    Ymodel.Tutsin_explicit_error.Data);

plot(Ymodel.Tutsin_implicit_error.Time,...
    Ymodel.Tutsin_implicit_error.Data);

plot(Ymodel.Halijak_error.Time,...
    Ymodel.Halijak_error.Data);

legend('Error between delayed Tutsin and continuous model',...
    'Error between implicit Tutsin and continuous model',...
    'Error between Halijak and continuous model');
grid minor;
xlabel('Time (s)');
ylabel('Output error');
title('Difference between discretized systems and continous model');

h_zod2=figure();
hold all
plot(Ymodel.Continuous.Time,Ymodel.Continuous.Data);
stairs(Ymodel.Tutsin_explicit.Time,...
    Ymodel.Tutsin_explicit.Data);

stairs(Ymodel.Tutsin_implicit.Time,...
    Ymodel.Tutsin_implicit.Data);

stairs(Ymodel.Halijak.Time,...
    Ymodel.Halijak.Data);

legend('Continuous model',...
    'Tutsin explicit',...
    'Tutsin implicit',...
    'Halijak');
grid minor;
xlabel('Time (s)');
ylabel('System response');
title('Outputs of continous system and discretized systems');

h_zod3=figure();
hold all
plot(Ymodel.Tutsin_implicit_effect.Time,...
    Ymodel.Tutsin_implicit_effect.Data);

legend('Difference between implicit and delayed Tutsin');
grid minor;
xlabel('Time (s)');
ylabel('Errors difference');
title('Difference between Tutsin output and Delayed Tutsin output');

h_zod4=figure();
hold all
plot(Ymodel.Continuous_command.Time,Ymodel.Continuous_command.Data);
plot(Ymodel.Tutsin_implicit_command.Time,...
    Ymodel.Tutsin_implicit_command.Data);

plot(Ymodel.Tutsin_explicit_command.Time,...
    Ymodel.Tutsin_explicit_command.Data);

plot(Ymodel.Halijak_command.Time,...
    Ymodel.Halijak_command.Data);

legend('Command from continuous model',...
    'Command from Tutsin',...
    'Command from delayed Tutsin',...
    'Command from Halijak');
grid minor;
xlabel('Time (s)');
ylabel('Command signals');
title('Commands of continous system and discretized systems');

%Saving figures
wPaperPos = [0 0 5 5];
wPaperSize = [5 5];

set(h_zod, 'PaperPosition', wPaperPos);
set(h_zod, 'PaperSize', wPaperSize);
saveas(h_zod,'Helico-Errors','pdf')

set(h_zod2, 'PaperPosition', wPaperPos);
set(h_zod2, 'PaperSize', wPaperSize);
saveas(h_zod2,'Helico-Step Response','pdf')

set(h_zod3, 'PaperPosition', wPaperPos);
set(h_zod3, 'PaperSize', wPaperSize);
saveas(h_zod3,...
    'Helico-Difference between tutsin implicit and explicit',...
    'pdf')

set(h_zod4, 'PaperPosition', wPaperPos);
set(h_zod4, 'PaperSize', wPaperSize);
saveas(h_zod4,'Helico-Commands','pdf')

%Saving model
saveas(get_param(model,'Handle'),'Helico-Model.pdf');
saveas(get_param(strcat(model,'/Tutsin implicit'),'Handle'),...
    'Helico-Model Tutsin.pdf');
saveas(get_param(strcat(model,'/Tutsin explicit'),'Handle'),...
    'Helico-Model Tutsin delayed.pdf');
close_system(model,false);


%Print the transfer functions on command windows for saving:
clc;
disp('Continuous transfer function')
disp('System:')
wSystem.mGetTf
disp('Compensator:')
wCompensator.mGetTf
disp('Closed Loop:')
wSystem.mGetClosedLoop(wCompensator.mGetTf,'continuous')

disp('Tutsin discretization')
disp('System:')
wSystem.mGetDiscreteTf('tutsin')
disp('Compensator:')
wCompensator.mGetDiscreteTf('tutsin')
disp('Closed Loop:')
wSystem.mGetClosedLoop(wCompensator.mGetTf,'tutsin')

disp('Tutsin discretization with delay')
disp('System:')
wSystem.mGetRetardedDiscreteTf('tutsin',1)
disp('Compensator:')
wCompensator.mGetRetardedDiscreteTf('tutsin',1)
disp('Closed Loop:')
wSystem.mGetClosedLoop(...
    wCompensator.mGetRetardedDiscreteTf('tutsin',1),...
    'tutsin')

disp('Halijak discretization')
disp('System:')
wSystem.mGetDiscreteTf('halijak')
disp('Compensator:')
wCompensator.mGetDiscreteTf('halijak')
disp('Closed Loop:')
wSystem.mGetClosedLoop(wCompensator.mGetTf,'halijak')
