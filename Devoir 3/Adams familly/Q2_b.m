
clear all; clc

T = 0.05;   % Période d'echantillonnage
tStop = 1; % Durée totale de la simulation
t=linspace(0,tStop,tStop*(1/T));
 
%Conditins initiaux
f1c(1) = 0;   f2c(1) = 0; f3c(1) = 1;
x1c(1) = 0;   x2c(1) = 0; x3c(1) = 0; 

f1c(2) = 0;   f2c(2) = 0; f3c(2) = 1;
x1c(2) = 0;   x2c(2) = 0; x3c(2) = 0; 

for n = 2:tStop/T-1

    x1p(n+1) = x1c(n) + 1/2*T*(3*f1c(n)-3*f1c(n-1));   
    x2p(n+1) = x2c(n) + 1/2*T*(3*f2c(n)-3*f2c(n-1));     
    x3p(n+1) = x3c(n) + 1/2*T*(3*f3c(n)-3*f3c(n-1));      
  
    f1p(n+1) = x2p(n+1);
    f2p(n+1) = x3p(n+1);
    f3p(n+1) = -200*x1p(n+1) - 30*x2p(n+1) - 11*x1p(n+1) + 1;
 
    x1c(n+1) = x1c(n) + (T/2)*(f1p(n+1)+ f1c(n));   
    x2c(n+1) = x2c(n) + (T/2)*(f2p(n+1)+ f2c(n));   
    x3c(n+1) = x3c(n) + (T/2)*(f3p(n+1)+ f3c(n));   
       
    f1c(n+1) = x2c(n+1);
    f2c(n+1) = x3c(n+1);
    f3c(n+1) = -200*x1c(n+1) - 30*x2c(n+1) - 11*x1c(n+1) + 1;
       
end

hold on
stairs(t,100*x2c,'r');
hold on
grid on;


