clear all
close all

t=0:.001:2*pi;
k=1;

z=(k'*exp(i*t))';
g=[1 2 4 6 8]'

for i = 1:length(g)
    rk2a(:,i)=sqrt(g(i)).*sqrt(4*z+g(i)-4)/2-g(i)/2;
    rk2b(:,i)=-sqrt(g(i)).*sqrt(4*z+g(i)-4)/2-g(i)/2;
end
plot([real(rk2a);real(rk2b)],[imag(rk2a);imag(rk2b)])
grid minor
legend('g=1', 'g=2', 'g=4', 'g=6', 'g=8')
title('RK2-stiff')
xlabel('real(\lambda T)')
ylabel('imag(\lambda T)')