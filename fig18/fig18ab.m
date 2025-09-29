% two solution-case2
clear all
clc;
%params
N=6; tau=1; 
k1=4*N-1;k2=4*N+2;
tstep = 14000;  
deltat = 1/tstep; 
t = 80;                                                                                                                           
tlist = 0:deltat:t; 
Cx = zeros(size(tlist)); Ce = zeros(size(tlist)); 
Cx(1) = sqrt(0.8); Ce(1) = sqrt(0.2); 
omegal=2.5*pi; Omega_3ls=2.8*pi; omegae=0.2*pi;
% Cx(1) = 1; Ce(1) = 0; 
% omegal=0;Omega_3ls=0;omegae=0;
omegax=Omega_3ls^2/(omegae+omegal-2*k1*pi/(N*tau))-cot(k1*pi/N)*(Omega_3ls^2/(omegae+omegal-2*k1*pi/(N*tau))-Omega_3ls^2/(omegae+omegal-2*k2*pi/(N*tau))+2*k1*pi/(N*tau)-2*k2*pi/(N*tau))/(cot(k1*pi/N)-cot(k2*pi/N))+2*k1*pi/(N*tau);
Gamma=2/N*(Omega_3ls^2/(omegae+omegal-2*k1*pi/(N*tau))-Omega_3ls^2/(omegae+omegal-2*k2*pi/(N*tau))+2*k1*pi/(N*tau)-2*k2*pi/(N*tau))/(cot(k1*pi/N)-cot(k2*pi/N));
abs((cot(k2*pi/N)*N*Gamma/2-Omega_3ls^2/(omegae+omegal-2*k2*pi/(N*tau))-omegal/2)/omegax)                                                     
sj=Gamma*t; 
sj1=fix(sj/10)*10;
s = N*Gamma/2 + 1i*(omegax-omegal/2);                      
for j=1:length(Cx)-1
    dCx1 = -s*Cx(j)-1i*Omega_3ls*Ce(j);                                                          
    dCe1=-1i*(omegae+omegal/2)*Ce(j)-1i*Omega_3ls*Cx(j); 
    
    dCx2 = -s*(Cx(j)+dCx1*deltat/2)-1i*Omega_3ls*(Ce(j)+dCe1*deltat/2);                                            
    dCe2=-1i*(omegae+omegal/2)*(Ce(j)+dCe1*deltat/2)-1i*Omega_3ls*(Cx(j)+dCx1*deltat/2); 
    
    dCx3 = -s*(Cx(j)+dCx2*deltat/2)-1i*Omega_3ls*(Ce(j)+dCe2*deltat/2);                                           
    dCe3=-1i*(omegae+omegal/2)*(Ce(j)+dCe2*deltat/2)-1i*Omega_3ls*(Cx(j)+dCx2*deltat/2); 
    
    dCx4 = -s*(Cx(j)+dCx3*deltat)-1i*Omega_3ls*(Ce(j)+dCe3*deltat/2);                                             
    dCe4=-1i*(omegae+omegal/2)*(Ce(j)+dCe3*deltat/2)-1i*Omega_3ls*(Cx(j)+dCx3*deltat/2); 
    DeltaCx = (dCx1+dCx2*2+dCx3*2+dCx4)/6;                                          
    DeltaCe = (dCe1+dCe2*2+dCe3*2+dCe4)/6; 
    
    for m=1:N
        for n=1:N                                                          
            index1 = j - abs(m-n)*tau*tstep;                               %t-|m-n|*tau
            HH1=exp(1i*omegal/2*abs(m-n));
            if tlist(j)-abs(m-n)*tau>0                                     %|m-n|*tau<t
                DeltaCx = DeltaCx - (Gamma)/2*HH1^tau*Cx(index1);
                
            elseif tlist(j)-(m+n)*tau>0                                    %t>(m+n)*tau
                DeltaCx = DeltaCx - (Gamma)/2*HH1^tau*Cx(index1)
            
            else
                DeltaCx = DeltaCx;
            end
        end     
    end
   DeltaCx = DeltaCx + N*(Gamma)/2*Cx(j);  
   Cx(j+1)=Cx(j)+DeltaCx*deltat;                                             %b(t+deltat)=b(t)+Deltab*deltat
   Ce(j+1)=Ce(j)+DeltaCe*deltat;
end
%% analytical  solution
Delta1=omegax-omegal/2;
Delta2=omegae+omegal/2;
ce=zeros(size(tlist));cx=zeros(size(tlist));
omegak1=2*k1*pi/(N*tau)-omegal/2;
omegak2=2*k2*pi/(N*tau)-omegal/2;
S_k1=-1i*omegak1;
S_k2=-1i*omegak2;
 for i=1:length(tlist)
    t=tlist(i);
    ce(i)=(Ce(1)*(S_k1+1i*Delta1+0.5*Gamma*(2*N/(1-exp(1i*2*k1*pi/N))-N))-1i*Omega_3ls*Cx(1))*exp(S_k1*t)...
        /((1+Gamma/2*(N*tau/((sin(k1*pi/N))^2)))*(S_k1+1i*Delta2)+(S_k1+1i*Delta1+Gamma/2*(2*N/(1-exp(1i*2*k1*pi/N))-N)))+(Ce(1)*(S_k2+1i*Delta1+0.5*Gamma*(2*N/(1-exp(1i*2*k2*pi/N))-N))-1i*Omega_3ls*Cx(1))*exp(S_k2*t)...
        /((1+Gamma/2*(N*tau/((sin(k2*pi/N))^2)))*(S_k2+1i*Delta2)+(S_k2+1i*Delta1+Gamma/2*(2*N/(1-exp(1i*2*k2*pi/N))-N)));
    
    cx(i)=(Cx(1)*(S_k1+1i*Delta2)-1i*Omega_3ls*Ce(1))*exp(S_k1*t)...
        /((1+Gamma/2*(N*tau/((sin(k1*pi/N))^2)))*(S_k1+1i*Delta2)+(S_k1+1i*Delta1+Gamma/2*(2*N/(1-exp(1i*2*k1*pi/N))-N)))+(Cx(1)*(S_k2+1i*Delta2)-1i*Omega_3ls*Ce(1))*exp(S_k2*t)...
        /((1+Gamma/2*(N*tau/((sin(k2*pi/N))^2)))*(S_k2+1i*Delta2)+(S_k2+1i*Delta1+Gamma/2*(2*N/(1-exp(1i*2*k2*pi/N))-N)));
 end   
Tlist=Gamma*tlist;
Y1=[abs(Cx).^2;abs(cx).^2];
% 
subplot('position',[0.12 0.61 0.366471449487555 0.3788461538461154])

plot1 = plot(Tlist,Y1);
set(plot1(1),'DisplayName','Eq. (43)','LineWidth',1.1,'Color',[64/255 123/255 208/255]);
set(plot1(2),'DisplayName','Eq. (62)','MarkerSize',1.3,'LineWidth',2.3,...
    'LineStyle','-.',...
    'Color',[163/255 42/255 49/255]);

xlim([0,40])
%xticks(linspace(0,sj1,6))
%xticklabels(linspace(0,,6))
xlabel('$\Gamma t$','Interpreter','latex')
ylabel('$|C_x(t)|^{2}$','Interpreter','latex')
%ylim([0 1.05]);

% 
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
%pbaspect([4 2.5 1])

% 
lgd = legend;
lgd.FontSize = 10;
%lgd.NumColumns = 2;
lgd.Interpreter = 'latex';
legend('boxoff');
title('(a)','position',[2.5,0.55],'FontSize',14);
hold on
% 
% axes2 = axes('Position',[0.09 0.25 0.4 0.55]);
% hold(axes2,'on');
% plot1 = plot(Tlist,Y1);
% set(plot1(1),'LineWidth',1.1,'Color',[64/255 123/255 208/255]);
% set(plot1(2),'LineWidth',1.3,...
%     'LineStyle','-.',...
%     'Color',[163/255 42/255 49/255]);
% % 
%  xlim(axes2,[36 40]);
% % 
%  ylim(axes2,[0 0.45]);
% box(axes2,'on');
% hold(axes2,'off');
% 
% set(axes2,'FontName','Times New Roman','FontSize',10,'LineWidth',1,...
%     'PlotBoxAspectRatio',[2 1 1]);
% hold(axes2,'on');
%%%%Ce   
% 
Y2=[abs(Ce).^2;abs(ce).^2];
subplot('position',[0.60 0.61 0.366471449487555 0.3788461538461154])

plot1 = plot(Tlist,Y2);
set(plot1(1),'DisplayName','Eq. (44)','LineWidth',1.1,'Color',[64/255 123/255 208/255]);
set(plot1(2),'DisplayName','Eq. (63)','MarkerSize',1.3,'LineWidth',2.3,...
    'LineStyle','-.',...
    'Color',[163/255 42/255 49/255]);

xlim([0,40])
xlabel('$\Gamma t$','Interpreter','latex')
ylabel('$|C_e(t)|^{2}$','Interpreter','latex')

% 
 %xlim(axes1,[0 50]);
% 
% ylim([0 1.05]);

% 
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
%pbaspect([4 2.5 1])
% 
lgd = legend;
%lgd.NumColumns = 2;
lgd.FontSize = 10;
lgd.Interpreter = 'latex';
legend('boxoff');
title('(b)','position',[2.5,0.13],'FontSize',14);
hold on
% 
% axes2 = axes('Position',[0.327994791666665 0.31972251867663 0.0958999999999986 0.1364]);
% hold(axes2,'on');
% plot1 = plot(Tlist,Y2);
% set(plot1(1),'LineWidth',1.1,'Color',[79/255 79/255 79/255]);
% set(plot1(2),'LineWidth',1.3,...
%     'LineStyle','-.',...
%     'Color',[152/255 245/255 255/255]);
% % 
%  xlim(axes2,[36 40]);
% % 
%  ylim(axes2,[0 0.1]);
% box(axes2,'on');
% hold(axes2,'off');
% % 
% set(axes2,'FontName','Times New Roman','FontSize',10,'LineWidth',1,...
%     'PlotBoxAspectRatio',[2 1 1]);
% hold(axes2,'on');

