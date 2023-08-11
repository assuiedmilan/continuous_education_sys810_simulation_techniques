addpath(genpath('..\..'));
dbstop if error

% ************************************************************ %
% ****************      INITIALIZATION      ****************** %
% ************************************************************ %

close all;
clear all; %#ok<CLSCR>
clc;

wSampleTimes=[0.3,0.2,0.19,0.18,0.17,0.16,0.15,0.1];
wSimulationTime=10;

wAdamsBashforth = tf([3,-1],[2,-2,0]);
wModifiedTutsin = tf([1 1 0],[2 -2 0]);

wContinuousSystemNum = [100,0];
wContinuousSystemDen = [1,11,30,200];

wPloter = Ploter([0 0 5 5],[5 5]);

for i=1:length(wSampleTimes)
    
    close all;
    java.lang.Runtime.getRuntime.freeMemory;
    
    clearvars -except wSampleTimes wSystem wSimulationTime wContinuousSystemNum...
        wContinuousSystemDen wAdamsBashforthNum wAdamsBashforthDen wPloter...
        wAdamsBashforth wModifiedTutsin i eval(SimOutput)
    
    wSampleTime = wSampleTimes(i);
    wMaxStep=wSampleTime/1000;
    
    wSystem = Discretizer(wSampleTime,...
        wContinuousSystemNum,...
        wContinuousSystemDen);
    
    adamsFamillyModelSetup;
    
    try
        load(wSaveFileName);
        Y = eval(wSaveFileName);
    catch
        h = errordlg(['Could not load', wSaveFileName, '. Simulation to be executed once. Exiting.']);
        waitfor(h);
        break;
    end
    
    % ************************************************************ %
    % ****************  PREDICTION - CORRECTION ****************** %
    % ************************************************************ %
    
    %Parameters
    clear t f1c f2c f3c x1c x2c x3c f1p f2p f3p x1p x2p x3p;
    t=linspace(0,wSimulationTime,wSimulationTime*(1/wSampleTime));
    
    %Stability region
    wSystem.mComputeStabilityRegion('Tutsin-Adam-Brashforth',wAdamsBashforth,wModifiedTutsin);
    %Initial conditions
    f1c(1) = 0;   f2c(1) = 0; f3c(1) = 1;
    x1c(1) = 0;   x2c(1) = 0; x3c(1) = 0;
    
    for n = 0:wSimulationTime/wSampleTime-2
        
        if (n == 0)
            x1p(n+2) = x1c(n+1) + (wSampleTime/2)*(3*f1c(n+1)); %#ok<*SAGROW>
            x2p(n+2) = x2c(n+1) + (wSampleTime/2)*(3*f2c(n+1));
            x3p(n+2) = x3c(n+1) + (wSampleTime/2)*(3*f3c(n+1));
        else
            x1p(n+2) = x1c(n+1) + (wSampleTime/2)*(3*f1c(n+1)-f1c(n));
            x2p(n+2) = x2c(n+1) + (wSampleTime/2)*(3*f2c(n+1)-f2c(n));
            x3p(n+2) = x3c(n+1) + (wSampleTime/2)*(3*f3c(n+1)-f3c(n));
        end
        
        f1p(n+2) = x2p(n+2);
        f2p(n+2) = x3p(n+2);
        f3p(n+2) = -200*x1p(n+2) - 30*x2p(n+2) -11*x3p(n+2) + 1;
        
        x1c(n+2) = x1c(n+1) + (wSampleTime/2)*(f1p(n+2)+ f1c(n+1));
        x2c(n+2) = x2c(n+1) + (wSampleTime/2)*(f2p(n+2)+ f2c(n+1));
        x3c(n+2) = x3c(n+1) + (wSampleTime/2)*(f3p(n+2)+ f3c(n+1));
        
        f1c(n+2) = x2c(n+2);
        f2c(n+2) = x3c(n+2);
        f3c(n+2) = -200*x1c(n+2) - 30*x2c(n+2) -11*x3c(n+2) + 1;
        
    end
    
    % ************************************************************ %
    % *******************  EXTRACTING FIGURES ******************** %
    % ************************************************************ %
    
    %Drawing and saving Plots
    
    wPloter.mDrawStandardPlot({[Y.Observable_continuous.Time...
        ,Y.Observable_continuous.Data]...
        ,[t;100*x2c]}...
        ,'stairs'...
        ,['Prediction-Correction vs Continuous ',strrep(num2str(wSampleTime*1000),'.',''),'ms']...
        ,'Time (s)'...
        ,'Step Response'...
        ,{'Commandable continuous';'Prediction-Correction'});
    
    close all;
end