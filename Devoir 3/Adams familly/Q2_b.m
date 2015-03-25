
clear all , clc, close all

T = 0.01;   % Période d'echantillonnage
tStop = 10; % Durée totale de la simulation
t=linspace(0,tStop,tStop*(1/T));

%Conditins initiaux
f1c(1) = 0;   f2c(1) = 0; f3c(1) = 1;
x1c(1) = 0;   x2c(1) = 0; x3c(1) = 0;

for n = 0:tStop/T-2
    
    if (n == 0)        
        x1p(n+2) = x1c(n+1) + (T/2)*(f1c(n+1));
        x2p(n+2) = x2c(n+1) + (T/2)*(f2c(n+1));
        x3p(n+2) = x3c(n+1) + (T/2)*(f3c(n+1));
    else        
        x1p(n+2) = x1c(n+1) + (T/2)*(f1c(n+1)-f1c(n));
        x2p(n+2) = x2c(n+1) + (T/2)*(f2c(n+1)-f2c(n));
        x3p(n+2) = x3c(n+1) + (T/2)*(f3c(n+1)-f3c(n));
    end
    
    f1p(n+2) = x2p(n+2);
    f2p(n+2) = x3p(n+2);
    f3p(n+2) = -200*x1p(n+2) - 30*x2p(n+2) -11*x3p(n+2) + 1;
    
    x1c(n+2) = x1c(n+1) + (T/2)*(f1p(n+2)+ f1c(n+1));
    x2c(n+2) = x2c(n+1) + (T/2)*(f2p(n+2)+ f2c(n+1));
    x3c(n+2) = x3c(n+1) + (T/2)*(f3p(n+2)+ f3c(n+1));
    
    f1c(n+2) = x2c(n+2);
    f2c(n+2) = x3c(n+2);
    f3c(n+2) = -200*x1c(n+2) - 30*x2c(n+2) -11*x3c(n+2) + 1;
    
end

figure()
hold on
stairs(t,100*x2c,'r');
hold on
grid on;


