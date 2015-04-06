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

wSampleTime = 0.001;
wSimulationTime = 5;

wSystem = Discretizer(wSampleTime,...
    A,B,C,D);

% ************************************************************ %
% *****************      RK3 STABILITY      ****************** %
% ************************************************************ %

wMatlabIndexBias = 1;
wRungeKuttaOrder = 4;
wPonderingCoefficients = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05];
wStabilityPolynom = zeros(1,wRungeKuttaOrder+1);

wZvalues = (exp(1i*(0:0.01:2*pi)))';
wStabilityValues = zeros(wRungeKuttaOrder*length(wZvalues),1);
wPlotData = zeros(length(wZvalues)*wRungeKuttaOrder*length(wPonderingCoefficients),2);

for k=1:wRungeKuttaOrder
    wStabilityPolynom(k+wMatlabIndexBias) = 1/factorial(k);
end
%Correctly order coefficients for matlab
wStabilityPolynom = fliplr(wStabilityPolynom);

for l = 1:length(wPonderingCoefficients)
    for k = 1:length(wZvalues)

        wStabilityPolynom(wRungeKuttaOrder+1) = 1-wPonderingCoefficients(l)*wZvalues(k);
        wRoots = roots(wStabilityPolynom);

        for h=0:wRungeKuttaOrder-1
            wStabilityValues(wRungeKuttaOrder*k-h,l) = wRoots(h+wMatlabIndexBias);
        end

    end
    
    wPlotData((l-1)*length(wZvalues)*wRungeKuttaOrder + wMatlabIndexBias:l*length(wZvalues)*wRungeKuttaOrder,1:2) = [real(wStabilityValues(:,l)),imag(wStabilityValues(:,l))];
end


 wPloter.mDrawStandardPlot({{wPlotData,'.'}}...
    ,'plot'...
    ,'Stability Region Runge-Kutta order 4'...
    ,'Real axis'...
    ,'Imaginary axis'...
    ,'Stability region');



% ************************************************************ %
% *************      RK3 STIFFED STABILITY      ************** %
% ************************************************************ %

clear wStabilityPolynom wStabilityValues wPlotData
wStabilityPolynom = zeros(1,wRungeKuttaOrder+1);
wStabilityValues = zeros(wRungeKuttaOrder*length(wZvalues),1);
wPlotData = zeros(length(wZvalues)*wRungeKuttaOrder*length(wPonderingCoefficients),2);

g = 12;

for k=1:wRungeKuttaOrder
    if (k==wRungeKuttaOrder)
        wStabilityPolynom(k+wMatlabIndexBias) = 1/g;
    else
        wStabilityPolynom(k+wMatlabIndexBias) = 1/factorial(k);
    end
end
%Correctly order coefficients for matlab
wStabilityPolynom = fliplr(wStabilityPolynom);

for l = 1:length(wPonderingCoefficients)
    for k = 1:length(wZvalues)

        wStabilityPolynom(wRungeKuttaOrder+1) = 1-wPonderingCoefficients(l)*wZvalues(k);
        wRoots = roots(wStabilityPolynom);

        for h=0:wRungeKuttaOrder-1
            wStabilityValues(wRungeKuttaOrder*k-h,l) = wRoots(h+wMatlabIndexBias);
        end

    end
    
    wPlotData((l-1)*length(wZvalues)*wRungeKuttaOrder + wMatlabIndexBias:l*length(wZvalues)*wRungeKuttaOrder,1:2) = [real(wStabilityValues(:,l)),imag(wStabilityValues(:,l))];
end


wPloter.mDrawStandardPlot({{wPlotData,'.'}}...
    ,'plot'...
    ,['Stability Region Runge-Kutta order 4 adapted for stiffed systems with g=',num2str(g)]...
    ,'Real axis'...
    ,'Imaginary axis'...
    ,'Stability region');
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
% ********************      RK3 SRP      ********************* %
% ************************************************************ %

% ************************************************************ %
% *****************      RK3 EQUATIONS      ****************** %
% ************************************************************ %
t=linspace(0,wSimulationTime,wSimulationTime*(1/wSampleTime));

wNmax = wSimulationTime/wSampleTime-1;
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