addpath(genpath('..\..'));
dbstop if error

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





