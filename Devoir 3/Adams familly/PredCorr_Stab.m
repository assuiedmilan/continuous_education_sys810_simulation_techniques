addpath(genpath('..\..'));

clear all
close all
clc;

wAdamsBashforthNum = [3,-1];
wAdamsBashforthDen = [2,-2,0];

wTutsinNum = [1 1 0];
wTutsinDen = [2 -2 0];

wResult=Discretizer.mProcessStabilityRoots_PredictorCorrector(tf(wAdamsBashforthNum,wAdamsBashforthDen),tf(wTutsinNum,wTutsinDen));

wt = 0:0.0001:2*pi;
z=exp(i*wt);

wRoot1 = zeros(size(z));
wRoot2 = zeros(size(z));

for i=1:length(z)
    
wRoot1(i) = -(2*(z(i) + (z(i)*(3*z(i)^2 - 3*z(i) + 1))^(1/2)))/(3*z(i) - 1);
wRoot2(i) = -(2*(z(i) - (z(i)*(3*z(i)^2 - 3*z(i) + 1))^(1/2)))/(3*z(i) - 1);

end

figure();
plot([real(wRoot1) real(wRoot2)],[imag(wRoot1) imag(wRoot2)]);

figure();
plot(real(wRoot1),imag(wRoot1));
figure();
plot(real(wRoot2),imag(wRoot2));

addpath(genpath('..\..'));

wAdamsBashforthNum = 1;
wAdamsBashforthDen = [1 -1];

wTutsinNum = [1 0];
wTutsinDen = [1 -1];

wResult=Discretizer.mProcessStabilityRoots_PredictorCorrector(tf(wAdamsBashforthNum,wAdamsBashforthDen),tf(wTutsinNum,wTutsinDen));

wRoot1 = zeros(size(z));
wRoot2 = zeros(size(z));

for i=1:length(z)
    
wRoot1(i) = - (4*z(i) - 3)^(1/2)/2 - 1/2;
wRoot2(i) =  (4*z(i) - 3)^(1/2)/2 - 1/2;

end

figure();
plot([real(wRoot1) real(wRoot2)],[imag(wRoot1) imag(wRoot2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R�gion de stabilit�
FE_Tustin1 = (sqrt(2*z-1) - 1);  % pole -1+j
FE_Tustin2 = (-sqrt(2*z-1) - 1); % pole -1-j
figure();
plot(real(FE_Tustin1),imag(FE_Tustin1))
hold on
plot(real(FE_Tustin2),imag(FE_Tustin2))
grid on
title('Pr�dicteur correcteur')


wAdamsBashforthNum = 1;
wAdamsBashforthDen = [1 -1];

wTutsinNum = [1 1];
wTutsinDen = [2 -2];

wResult=Discretizer.mProcessStabilityRoots_PredictorCorrector(tf(wAdamsBashforthNum,wAdamsBashforthDen),tf(wTutsinNum,wTutsinDen));

wRoot1 = zeros(size(z));
wRoot2 = zeros(size(z));

for i=1:length(z)
    
wRoot1(i) = - (2*z(i) - 1)^(1/2) - 1;
wRoot2(i) =  (2*z(i) - 1)^(1/2) - 1;

end

figure();
plot([real(wRoot1) real(wRoot2)],[imag(wRoot1) imag(wRoot2)]);