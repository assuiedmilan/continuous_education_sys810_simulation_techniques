clear all
close all

z=tf('z')
wt=0:.01:2*pi;
k=[1];

z=(k'*exp(i*wt))';

subplot(2,1,1)
FE=(z-1);
plot(real(FE),imag(FE))
axis([-2 2 -1 1])
grid on
title('Euler avant')

subplot(2,1,2)
BE=(z-1)./z;
plot(real(BE),imag(BE))
axis([-2 2 -1 1])
grid on
title('Euler arri?re')

figure
subplot(2,1,1)
TUS=2*(z-1)./(z+1);
plot(real(TUS),imag(TUS))
axis([-1 1 -1 1])
grid on
title('Tustin')

subplot(2,1,2)
ab2=(2*z.*(z-1))./(3*z-1);
plot(real(ab2),imag(ab2))
axis([-1 1 -1 1])
grid on
title('AB2')

figure
subplot(2,1,1)
ab3=(12*z.^2.*(z-1))./(23*z.^2-16*z+5);
plot(real(ab3),imag(ab3))
axis([-1 1 -1 1])
grid on
title('AB3')

subplot(2,1,2)
ab4=(24*z.^3.*(z-1))./(55*z.^3-59*z.^2+37*z-9);
plot(real(ab4),imag(ab4))
axis([-1 1 -1 1])
grid on
title('AB4')