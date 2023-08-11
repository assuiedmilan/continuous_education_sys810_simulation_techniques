k2s=[1:1:50];
bs=[1:1:50];

%Specific 3-D post treatment

K2=zeros(length(k2s),1);
B=zeros(length(bs),1);
Xm=zeros(length(k2s),length(bs));

for i=1:length(k2s)
    k2 = k2s(i);
    K2(i)=k2;
    
    for j=1:length(bs)
        b=bs(j);
        B(j)=b;
        
        wSaveFileName     = strcat('masseRessort_k2_',..
			strrep(strcat(num2str(k2),'_b_',num2str(b)),'.','dec'));
        wSwitchOutputName = strcat('masseRessort_Switch_',wSaveFileName);
        load(wSaveFileName);
        wStruct = eval(wSwitchOutputName);
        
        Xm(i,j)=x0-min(wStruct.x_switch.Data);
        
         if (k2==1 && b==1)
             
             wIndexMin = find(wStruct.x_switch.Data==min(wStruct.x_switch.Data));
             
             figure();
             plot(wStruct.x_switch.Time,wStruct.x_switch.Data);
             hold('on');
             
             plot(wStruct.dx_switch.Time,wStruct.dx_switch.Data);
             plot(wStruct.Fs_switch.Time,wStruct.Fs_switch.Data);
             plot(wStruct.x_switch.Time(wIndexMin),..
				wStruct.x_switch.Data(wIndexMin),'ro');
             grid minor;
             
             title(strcat('Deplacement de la masse et efforts',..
				'appliques au cours du temps');
				(k_2=1 et b=1)');
             h = legend('x(t)','$\frac{dx}{dt}$','F_s(t)');
             set(h,'interpreter','latex');
             set(h,'position',[0.8 0.8 0.15 0.15]);
             set(gca,'fontsize',15)
             xlabel('Time (s)');
             
         end
        
    end
end

figure();
surf(K2,B,Xm);
grid minor;
legend('x_{Min}(k2,b)');
title('Deplacement maximal de la masse selon l''amortissement');
xlabel('K_2');
ylabel('b');
zlabel('Deviation maximale');

%Calcul du gradient
Gxm = gradient(Xm,1,1);

figure();
surf(K2,B,Xm);
grid minor;
h = legend('$\vec\nabla(X_{Min}(k2,b))$');
set(h,'interpreter','latex');
title('Gradient du deplacement maximal de la masse selon l''amortissement');
xlabel('K_2');
ylabel('b');
zlabel('Gradient');

wTreshold = 0.01
[r,c]=find(abs(Gxm)<wTreshold);
A=[r,c];

figure();
plot(r,c);
h = legend(strcat('Lieu des gradients < ',num2str(wTreshold)));
xlabel('K_2');
ylabel('b');