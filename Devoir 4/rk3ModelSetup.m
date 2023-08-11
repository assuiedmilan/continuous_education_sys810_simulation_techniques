addpath(genpath('..\..'));
dbstop if error

% ************************************************************ %
% ****************  SIMULINK INITIALIZATION ****************** %
% ************************************************************ %

model='rk3Model';
load_system(model)
tic

wSaveFileName     = 'SimOutput';
wContBlock = strcat(model,'/State-Space');

set_param(model,'StopFcn','save(wSaveFileName,wSaveFileName)');
set_param(strcat(model,'/Output'),'VariableName',wSaveFileName);


%Parametrage simulation discrète
set_param(model,'Solver','Ode3');
set_param(model,'FixedStep','wSampleTime');
set_param(wContBlock,'A','A','B','B','C','C','D','D','X0','X0');

set_param(model, 'StopTime', 'wSimulationTime');

set_param(model, 'MaxStep', 'wMaxStep');
