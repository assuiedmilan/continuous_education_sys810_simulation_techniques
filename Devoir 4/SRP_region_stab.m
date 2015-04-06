clear all
close all

t=0:.01:2*pi;
k=1;
z=(k'*exp(i*t))';

subplot(3,1,1)
a0 = 0.9; a1=-1.9; b0=-.09; b1=0.1899959;
srp110 = (z.^2+a1*z+a0)./(b1*z+b0);
plot(real(srp110),imag(srp110),'b')
grid on
title('srp110')
xlabel('real(\lambda T)')
ylabel('imag(\lambda T)')

subplot(3,1,2)
a0 = 0.99; a1=-1.99; b0=-.0099; b1=0.01899959;
srp1100 = (z.^2+a1*z+a0)./(b1*z+b0);
plot(real(srp1100),imag(srp1100),'r')
grid on
title('srp1100')
xlabel('real(\lambda T)')
ylabel('imag(\lambda T)')

subplot(3,1,3)
a0 = 0.999; a1=-1.999; b0=-.000999; b1=0.001899959;
srp11000 = (z.^2+a1*z+a0)./(b1*z+b0);
plot(real(srp11000),imag(srp11000),'k')
grid on
title('srp11000')
xlabel('real(\lambda T)')
ylabel('imag(\lambda T)')