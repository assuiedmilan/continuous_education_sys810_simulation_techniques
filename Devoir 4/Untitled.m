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

wSampleTimes = [0.05,0.04,0.03,0.02,0.015,0.01,0.005];
wStiffAdapters = [6,9,12,15];
wIdealSimulationTime = 5;
wRungeKuttaOrder = 3;
wSystem = Discretizer(1,A,B,C,D);

% ************************************************************ %
% *****************      RK3 STIFFNESS      ****************** %
% ************************************************************ %
wStiffnessRatio = max(abs(real(wSystem.mGetPoles('continuous'))))/min(abs(real(wSystem.mGetPoles('continuous'))));

%Stabilit� � partir de 20ms
%Pr�cision � partir de 10ms
for wSampleTimeIndex = 1:length(wSampleTimes)
    
    close all;
    java.lang.Runtime.getRuntime.freeMemory;
    
    clearvars -except wSampleTimes wStiffAdapters wIdealSimulationTime...
        wRungeKuttaOrder wSystem wPloter wSampleTimeIndex A B C D X0
        
    wSampleTime = wSampleTimes(wSampleTimeIndex)
    wSimulationTime = ceil(wIdealSimulationTime/wSampleTime)*wSampleTime;
    wSystem.mSetSampleTime(wSampleTime);    
    
    wSampleTimeLegend = [' ',strrep(num2str(wSampleTime*1000),'.',''),'ms'];
    
    % ************************************************************ %
    % **************** CONTINUOUS  SIMULATION ******************** %
    % ************************************************************ %
    rk3ModelSetup;
    
    myopts=simset('SrcWorkspace','current','DstWorkspace','current');
    
    sim(model,wSimulationTime,myopts);
    while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
    end
    
    t_sim = toc;
    fprintf('\nTemps de simulation => %3.3g s\n',t_sim)

    % ************************************************************ %
    % *****************      RK3 STABILITY      ****************** %
    % ************************************************************ %
    wSystem.mComputeStabilityRegion(['Runge-Kutta order ',num2str(wRungeKuttaOrder),wSampleTimeLegend],wRungeKuttaOrder);
    
    % ************************************************************ %
    % *************      RK3 STIFFED STABILITY      ************** %
    % ************************************************************ %
    wSystem.mComputeStabilityRegion(['Runge-Kutta order ',num2str(wRungeKuttaOrder),' precision ',num2str(wRungeKuttaOrder-1),wSampleTimeLegend],wRungeKuttaOrder,wStiffAdapters);
    
    % ************************************************************ %
    % *****************     RK3 DISCRETE POLES  ****************** %
    % ************************************************************ %    
    zExact  = zeros(1,length(wSampleTimes));
    zRk3    = zeros(1,length(wSampleTimes));
    zRk3p2  = zeros(1,length(wSampleTimes));
    
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
        ,['RK3 poles',wSampleTimeLegend]...
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
        ,['ODE3 vs RK3 ',wSampleTimeLegend]...
        ,'Time (s)'...
        ,'Step Response'...
        ,{'ODE3';'RK3 Recursive equations'});

    close all;
end
% ************************************************************ %
% *****************      SRP EQUATIONS      ****************** %
% ************************************************************ %
% wL = max(abs(wSystem.mGetPoles('continuous')));
% wT = wSystem.mGetSampleTime();
%
% E = [1 1 1 0 0 0;...
%     0 1 2 -1 -1 -1;...
%     0 1/2 2 0 -1 -2;...
%     1 0 0 -wL*wT 0 0;...
%     0 1 0 0 -wL*wT 0;...
%     0 0 1 0 0 -wL*wT];
%
% C = [-1;-3;-9/2;0;0;-exp(wL*wT)];
%
% X = linsolve(E,C);
%
% dbstop if error
% wIntegrator = tf([0,X(length(X)/2+1:length(X))'],[1,X(1:length(X)/2)']);
% wSystem.mComputeStabilityRegion('SRP',wIntegrator);
%
% t=0:.01:2*pi;
% k=1;
% z=(k'*exp(i*t))';
%
% srp110 = (z.^3+X(3)*z.^2+X(2)*z+X(1))./(X(6)*z.^2+X(5)*z+X(4));
% figure()
% plot(real(srp110),imag(srp110),'b')
%
