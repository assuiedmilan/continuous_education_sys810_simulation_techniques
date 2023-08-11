
wAi = cellstr(['Ai_0dec1';'Ai_1    ';'Ai_10   ']);
wAin = cellstr(['0.1';'1  ';'10 ']);

for i=1:size(wAi,1)
    
    wCurrentAi  = wAi(i);
    wCurrentAin = wAin(i);
    
    wCurrentAi  = wCurrentAi{1};
    wCurrentAin = wCurrentAin{1};
    
    figure();
    wFileName = strcat('circuit_',wCurrentAi,'.mat');
    load(wFileName);
    wLinearModelIl    =  eval(strcat('linear_circuit_',..
		wCurrentAi,'.linear_Il'));
    wNonlinearModelIl =  eval(strcat('nonLinear_circuit_',..
		wCurrentAi,'.NonLinear_model.nonLinear_Il'));
    wSimPowerIl       =  eval(strcat('nonLinear_circuit_',..
		wCurrentAi,'.SimPwr_model.simPwr_Il'));

    wTau = 0.67*mean(wLinearModelIl.Data);
    wIndexTau = find(wNonlinearModelIl.Data>wTau,1,'first');
    
    plot(wLinearModelIl.Time,wLinearModelIl.Data);
    hold('on');
    
    plot(wNonlinearModelIl.Time,wNonlinearModelIl.Data);
    plot(wSimPowerIl.Time,wSimPowerIl.Data);
    plot(wNonlinearModelIl.Time(wIndexTau),..
		wNonlinearModelIl.Data(wIndexTau),'ro');
    
    grid minor;
    legend(strcat('A=',wCurrentAin,'V, modele linearise'),..
		strcat('A=',wCurrentAin,'V, modele non lineaire'),..
		strcat('A=',wCurrentAin,'V, simPower'));
    xlabel('Time (s)');
    ylabel('Intensite (A)');
    
end

