addpath(genpath('..\..'));
dbstop if error

% ************************************************************ %
% ****************      INITIALIZATION      ****************** %
% ************************************************************ %

close all;
temporaryBreakpointData=dbstatus('-completenames');
clear functions; %#ok<CLSCR>
dbstop(temporaryBreakpointData);
clear global;
clear variables;
clc;

wPloter = Ploter([0 0 5 5],[5 5]);

A = [-112.5 53.3 42.6...
    ; -48.11 0.2050 47.92...
    ; 2.841 -1.269 -4.117];

B = [0.2733;0.4667;0.2667];
C = [1,0,0];
D = 0;
X0 = [-0.5; 0.7; -0.5];

wSampleTime = 0.020;
wSimulationTime = 10;

%Stabilité à partir de 20ms
%Précision à partir de 10ms
wSystem = Discretizer(wSampleTime,...
    A,B,C,D);

wRungeKuttaOrder = 3;
% ************************************************************ %
% *****************      RK3 STABILITY      ****************** %
% ************************************************************ %
wSystem.mComputeStabilityRegion(['Runge-Kutta order ',num2str(wRungeKuttaOrder)],wRungeKuttaOrder);

% ************************************************************ %
% *************      RK3 STIFFED STABILITY      ************** %
% ************************************************************ %
wStiffAdapter = [6,12];
wSystem.mComputeStabilityRegion(['Runge-Kutta order ',num2str(wRungeKuttaOrder),' with stiff adapter =',num2str(wStiffAdapter)],wRungeKuttaOrder,wStiffAdapter);

% ************************************************************ %
% **********************  SIMULATION ************************* %
% ************************************************************ %
rk3ModelSetup;

myopts=simset('SrcWorkspace','current','DstWorkspace','current');

sim(model,wSimulationTime,myopts);
while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
end

t_sim = toc;
fprintf('\nTemps de simulation => %3.3g s\n',t_sim)

% ************************************************************ %
% *****************      RK3 STIFFNESS      ****************** %
% ************************************************************ %
wStiffnessRatio = max(abs(real(wSystem.mGetPoles('continuous'))))/min(abs(real(wSystem.mGetPoles('continuous'))));

% ************************************************************ %
% *****************     RK3 DISCRETE POLES  ****************** %
% ************************************************************ %
wSampleTimes = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.001]
zExact = zeros(1,length(wSampleTimes));
zRk3 = zeros(1,length(wSampleTimes));
zRk3p2 = zeros(1,length(wSampleTimes));

for k=1:length(wSampleTimes)
    
wSystem.mSetSampleTime(wSampleTimes(k));
wLs = abs(wSystem.mGetPoles('continuous'));
wTs = wSystem.mGetSampleTime();
    
    for i=1:length(wLs)
    wLT = wLs(i)*wTs;
    zExact(k)  = abs(exp(wLT));
    zRk3(k)    = abs(1+wLT+1/2*wLT^2+1/6*wLT^3);
    zRk3p2(k)  = abs(1+wLT+1/2*wLT^2+1/9*wLT^3);
    end    
    
    
end

wPloter.mDrawStandardPlot({{[wSampleTimes;zExact],'.'}...
    ,{[wSampleTimes;zRk3],'*'}...
    ,{[wSampleTimes;zRk3p2],'o'}}...
    ,'plot'...
    ,'RK3 poles'...
    ,'Sample time (s)'...
    ,'|z|'...
    ,{'Exact pole';'RK3 pole';'RK3 precision 2 pole'});

% ************************************************************ %
% *****************      RK3 EQUATIONS      ****************** %
% ************************************************************ %
t=linspace(0,wSimulationTime,wSimulationTime*(1/wSampleTime));
wMatlabIndexBias = 1;
wNmax = fix(wSimulationTime/wSampleTime)-1;
U=ones(1,wNmax);
X = zeros(3,1);
X(:,1) = X0;
for n = 0:wNmax-1
    
    k1  = A*X(:,n+wMatlabIndexBias) + B*U(n+wMatlabIndexBias);
    t1  = 1/2 * wSampleTime;
    xp1 = X(:,n+wMatlabIndexBias) + t1 * k1;
    up1 = 1;
    
    k2  = A*xp1 + B*up1;
    t2  = 3/4 * wSampleTime;
    xp2 = X(:,n+wMatlabIndexBias) + t2 * k1;
    up2 = 1;
    
    k3 = A*xp2 + B*up2;
    
    X(:,n+wMatlabIndexBias+1) = X(:,n+wMatlabIndexBias) + 1/9*wSampleTime*(2*k1+3*k2+4*k3);
    
end

wPloter.mDrawStandardPlot({[SimOutput.Continuous_signal.Time...
    ,SimOutput.Continuous_signal.Data]...
    ,[t;C*X]}...
    ,'stairs'...
    ,['ODE3 vs RK3 ',strrep(num2str(wSampleTime*1000),'.',''),'ms']...
    ,'Time (s)'...
    ,'Step Response'...
    ,{'ODE3';'RK3 Recursive equations'});