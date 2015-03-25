
clear all , clc

T = 0.1;   % Période d'echantillonnage
tStop = 10; % Durée totale de la simulation
t=linspace(0,tStop,tStop*(1/T));
 
%Conditins initiaux
f1c(1) = 0;   f2c(1) = 1;
x1c(1) = 0;   x2c(1) = 0; 

for n = 1:tStop/T-1

    x1p(n+1) = x1c(n) + T*f1c(n);   
    x2p(n+1) = x2c(n) + T*f2c(n);   
  
    f1p(n+1) = x2p(n+1);
    f2p(n+1) = -2*x1p(n+1) - 2*x2p(n+1) + 1;
 
    x1c(n+1) = x1c(n) + (T/2)*(f1p(n+1)+ f1c(n));   
    x2c(n+1) = x2c(n) + (T/2)*(f2p(n+1)+ f2c(n));   
       
    f1c(n+1) = x2c(n+1);
    f2c(n+1) = -2*x1c(n+1) - 2*x2c(n+1) + 1;
       
end

hold on
stairs(t,x1c,'r');
hold on
grid on;


