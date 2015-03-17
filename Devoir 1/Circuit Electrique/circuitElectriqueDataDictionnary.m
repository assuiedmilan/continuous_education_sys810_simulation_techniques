addpath(genpath('..\..'));

clear all; %#ok<CLSCR>
%Input tension amplitude
Ais=[0.1,1,10];
Ei=2;

%Varistor parameters
Io=1/8;
Vo=1;
alpha=3;

%Circuit constants
Ca=2;
L=3;
R=1;

%Point d'equilibre
Ei_mean = Ei;
u0 = Ei_mean;
x0 = [Ei_mean;1/8*Ei_mean^3];

%Representation d'etat linearisee autour du point d'equilibre:
A=[-3/(8*Ca)*Ei_mean^2-1/(R*Ca), 1/Ca;-1/L,0];
B=[1/(R*Ca);1/L];
C=eye(2);
D=[0;0];

%Parametres de simulation
t0=0;
tf=20;
max_step=0.1;
min_step=0.001;
tol=0.001;

%Jacobian computation
syms x1 x2 ei x1_dot x2_dot Rj Cj Lj eij
x1_dot = -1/(8*Cj)*x1^3 - 1/(Rj*Cj)*x1 + 1/Cj*x2 + 1/(Rj*Cj)*eij;
x2_dot = -1/Lj*x1 + 1/Lj*eij;
Aj = jacobian([x1_dot;x2_dot],[x1;x2])
Bj = jacobian([x1_dot;x2_dot],eij)

As = subs(Aj,Cj,Ca);
As = subs(As,Rj,R);
As = subs(As,Lj,L);
As = subs(As,x1,Ei_mean);
As = double(As);

Bs = subs(Bj,Cj,Ca);
Bs = subs(Bs,Rj,R);
Bs = subs(Bs,Lj,L);
Bs = double(Bs);

%Parametrage simulation
model='circuitElectrique';
load_system(model)
tic
for i=1:length(Ais)
    Ai = Ais(i);
    
    wSaveFileName        = strcat('circuit_Ai_',..
		strrep(num2str(Ai),'.','dec'));
    wLinearOutputName    = strcat('linear_',wSaveFileName);
    wNonLinearOutputName = strcat('nonLinear_',wSaveFileName);
    
    set_param(model,'StopFcn',..
		'save(wSaveFileName,wLinearOutputName,wNonLinearOutputName)');
    set_param(strcat(model,'/linearOutput'),..
		'VariableName',wLinearOutputName);
    set_param(strcat(model,'/nonLinearOutput'),..
		'VariableName',wNonLinearOutputName);    
    
    myopts=simset('SrcWorkspace','current','DstWorkspace','current');
    
    sim(model,tf,myopts);
    while (strcmp(get_param(model,'SimulationStatus'),'stopped')==0);
    end
end
t_sim = toc;
close_system(model,false)
fprintf('\nTemps de simulation => %3.3g s\n',t_sim)

