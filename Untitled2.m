wt=0:.0001:2*pi;
x=exp(i*wt);

tic
FE_1 = (sqrt(2*x-1) - 1); 
FE_2 = (-sqrt(2*x-1) - 1); 
toc

plot(real(FE_1),imag(FE_1))
hold on
plot(real(FE_2),imag(FE_2))

z=sym('z');
FEz_1 = (sqrt(2*z-1) - 1); 
FEz_2 = (-sqrt(2*z-1) - 1); 

tic
FE_1 = subs(FEz_1,x);
FE_2 = subs(FEz_2,x);
toc

tic
FE_1 = double(FE_1);
FE_2 = double(FE_2);
toc