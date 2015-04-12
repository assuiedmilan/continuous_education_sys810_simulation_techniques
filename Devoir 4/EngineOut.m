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

wSampleTimes = [0.05,0.04,0.03,0.02,0.015,0.01,0.005,0.001];
wStiffAdapters = [6,9,12,15];
wIdealSimulationTime = 5;
wRungeKuttaOrder = 3;
wSystem = Discretizer(1,A,B,C,D);

% ************************************************************ %
% *****************      RK3 STIFFNESS      ****************** %
% ************************************************************ %
wStiffnessRatio = max(abs(real(wSystem.mGetPoles('continuous'))))/min(abs(real(wSystem.mGetPoles('continuous'))));

% ************************************************************ %
% *****************     RK3 DISCRETE POLES  ****************** %
% ************************************************************ %
wSampleTimesVector = 0.01:0.01:0.4;
wL = min(abs(wSystem.mGetPoles('continuous')));

zExact  = zeros(1,length(wSampleTimesVector));
zRk3    = zeros(1,length(wSampleTimesVector));
zRk3p2  = zeros(1,length(wSampleTimesVector));

for g=1:length(wStiffAdapters)
    
    for k=1:length(wSampleTimesVector)
        
        wLT = wL*wSampleTimesVector(k);
        
        zExact(k)  = abs(exp(wLT));
        zRk3(k)    = abs(1+wLT+1/2*wLT^2+1/6*wLT^3);
        zRk3p2(k)  = abs(1+wLT+1/2*wLT^2+1/wStiffAdapters(g)*wLT^3);
        
    end
    
    wPloter.mDrawStandardPlot({{[wSampleTimesVector;zExact],'.'}...
        ,{[wSampleTimesVector;zRk3],'*'}...
        ,{[wSampleTimesVector;zRk3p2],'o'}}...
        ,'plot'...
        ,['RK3 pole vs fastest exact pole with g= ',num2str(wStiffAdapters(g))]...
        ,'Sample time (s)'...
        ,'|z|'...
        ,{'Exact pole';'RK3 pole';'RK3 precision 2 pole'});
    
end

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
    % *****************      RK3 EQUATIONS      ****************** %
    % ************************************************************ %
    
    for h=1:1
        g = wStiffAdapters(h); %g=9 => stable at 40ms, precise at 0.15 ?
        
        %Kutta's third-order method
        %0   |  0   0   0
        %1/2 |  1/2 0   0
        %3/4 |  0   3/4 0
        %------------------
        %    | 2/9  3/9 4/9
        wNmax = fix(wSimulationTime/wSampleTime)-1;
        [oX,oY,oT] = RK3.sProcessRungeKutta([1/6,2/3,1/6], [1/2,0;-1,2], X0, ones(1,wNmax), {A,B,C,D}, wNmax, wSampleTime);
        
        wHandle = wPloter.mDrawStandardPlot({[SimOutput.Continuous_signal.Time...
            ,SimOutput.Continuous_signal.Data]...
            ,[oT;oY]}...
            ,'stairs'...
            ,['ODE3 vs RK3 g=',num2str(g),wSampleTimeLegend]...
            ,'Time (s)'...
            ,'Step Response'...
            ,{'ODE3';'RK3 Recursive equations'});
        RK3.sAddButcherTable(wHandle,[1/6,2/3,1/6], [1/2,0;-1,2]);
        
        wPloter.mDrawStandardPlot({[oT;oY]}...
            ,'stairs'...
            ,['RK3 Huen g=',num2str(g),wSampleTimeLegend]...
            ,'Time (s)'...
            ,'Step Response'...
            ,{'RK3 Recursive equations'});
        
    end
    wPloter.mDrawStandardPlot({[SimOutput.Continuous_signal.Time...
        ,SimOutput.Continuous_signal.Data]}...
        ,'stairs'...
        ,['ODE3',wSampleTimeLegend]...
        ,'Time (s)'...
        ,'Step Response'...
        ,{'ODE3'});
    close all;
end

%Saving Model
wPloter.mProcessSaveModel(model);

