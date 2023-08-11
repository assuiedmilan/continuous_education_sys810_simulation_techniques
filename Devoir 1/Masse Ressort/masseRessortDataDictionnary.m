addpath(genpath('..\..'));

clear all; %#ok<CLSCR>

%Parametres
g=9.8; %#ok<*NASGU>
m=1;
k1=1;
k2s=1:1:50;
bs=1:1:50;
b2=0;
v0=-3;
x1=1;
x0=2;

%Parametres de simulation
t0=0;
tf=20;
max_step=0.1;
min_step=0.0001;
tol=0.0001;

%Parametrage simulation
model='masseRessort';
load_system(model)
tic
for i=1:length(k2s)
    k2 = k2s(i);
    
    for j=1:length(bs)
        b=bs(j);
        
        wSaveFileName     = strcat('masseRessort_k2_',..
			strrep(strcat(num2str(k2),'_b_',num2str(b)),'.','dec'));
        wSwitchOutputName = strcat('masseRessort_Switch_',wSaveFileName);
        wMatlabOutputName = strcat('masseRessort_Matlalb_',wSaveFileName);
        
        set_param(model,'StopFcn',..
		'save(wSaveFileName,wSwitchOutputName,wMatlabOutputName)');
        set_param(strcat(model,'/switchOutput'),..
		'VariableName',wSwitchOutputName);
        set_param(strcat(model,'/matlabOutput'),..
		'VariableName',wMatlabOutputName);
        
        myopts=simset('SrcWorkspace','current',..
			'DstWorkspace','current');
        
        sim(model,tf,myopts);
        while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
        end
    end
    
end
t_sim = toc;
close_system(model,false)
fprintf('\nTemps de simulation => %3.3g s\n',t_sim)

