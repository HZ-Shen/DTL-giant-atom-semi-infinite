%% 
%%%
clear all
clc;
N=3; tau=1; 
tstep = 8000;                                                             
deltat = 1/tstep;                                                         
t = 120;                                                                                                                           
tlist = 0:deltat:t;   
Cx = zeros(size(tlist)); Ce = zeros(size(tlist)); 
% omegax=2*pi/tau;                                                          
% Gamma=0.05*pi/tau;
%%%new params
%omegae=0;omegal=0; Omega_3ls=0;  
%Cx(1) = 1; Ce(1) = 0; 
Cx(1) = sqrt(0.8); Ce(1) = sqrt(0.2); 
omegax=2*pi/tau; 
k=4; omegal=1.2*pi;Omega_3ls=-0.6*pi;Gamma=0.3*pi/tau; 
omegae=2*Omega_3ls^2/(Gamma*N*cot(k*pi/N)-4*k*pi/(N*tau)+2*omegax)-omegal+2*k*pi/(N*tau);%omegae has expression                                                    
abs((cot(k*pi/N)*N*Gamma/2-Omega_3ls^2/(omegae+omegal-2*k*pi/(N*tau))-omegal/2)/omegax) 
delta_w=0*Gamma;                                                        
Gamma_e=0*Gamma;
R=1;
r=R+1i*sqrt(R*(1-R));
s = N*(Gamma+Gamma_e)/2 + 1i*(omegax-omegal/2+sqrt(2*delta_w)); 
Delta1=omegax-omegal/2; Delta2=omegae+omegal/2;
ce=zeros(size(tlist));cx=zeros(size(tlist));
Omegak=2*k*pi/(N*tau)-omegal/2;                                            
S_k=-1i*Omegak;
condition=Cx(1)*(S_k+1*i*Delta2)-1i*Omega_3ls*Ce(1);
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
                DeltaCx = DeltaCx - (Gamma)/2*HH1.^tau*Cx(index1);
         
            else
                DeltaCx = DeltaCx;
            end
        end     
    end
   DeltaCx = DeltaCx + N*(Gamma)/2*Cx(j);  
   Cx(j+1)=Cx(j)+DeltaCx*deltat;                                             %b(t+deltat)=b(t)+Deltab*deltat
   Ce(j+1)=Ce(j)+DeltaCe*deltat;
end
%%%analytical solution
for i=1:length(tlist)
    t=tlist(i);
    ce(i)=(Ce(1)*(S_k+1i*Delta1+0.5*Gamma*(2*N/(1-exp(1i*2*k*pi/N))-N))-1i*Omega_3ls*Cx(1))*exp(S_k*t)...
        /((1+Gamma/2*(N*tau/((sin(k*pi/N))^2)))*(S_k+1i*Delta2)+(S_k+1i*Delta1+Gamma/2*(2*N/(1-exp(1i*2*k*pi/N))-N)));
    
    cx(i)=(Cx(1)*(S_k+1i*Delta2)-1i*Omega_3ls*Ce(1))*exp(S_k*t)...
        /((1+Gamma/2*(N*tau/((sin(k*pi/N))^2)))*(S_k+1i*Delta2)+(S_k+1i*Delta1+Gamma/2*(2*N/(1-exp(1i*2*k*pi/N))-N)));
 end 
Tlist=Gamma*tlist;
Y1=[abs(Cx).^2;abs(cx).^2];
% 
subplot('position',[0.60 0.61 0.366471449487555 0.378846153846154])
plot1 = plot(Tlist,Y1);
set(plot1(1),'DisplayName','Eq. (43)','LineWidth',1.5,'Color',[64/255 123/255 208/255]);
set(plot1(2),'DisplayName','Eq. (60)','MarkerSize',2,'LineWidth',2.3,...
    'LineStyle','-.',...
    'Color',[163/255 42/255 49/255]);

xlim([0,100])
%xticks(linspace(0,sj1,6))
%xticklabels(linspace(0,,6))
xlabel('$\Gamma t$','Interpreter','latex')
ylabel('$|C_x(t)|^{2}$','Interpreter','latex')
%ylim([0 1.05]);

% 
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
pbaspect([4 2.5 1])

% 
lgd = legend;
lgd.FontSize = 11;
lgd.Interpreter = 'latex';
legend('boxoff');

title('(b)','position',[6,0.66],'FontSize',14);

%%%%Ce   
% 
Y2=[abs(Ce).^2;abs(ce).^2];
subplot('position',[0.60 0.16 0.366471449487555 0.378846153846154])
plot1 = plot(Tlist,Y2);
set(plot1(1),'DisplayName','Eq. (44)','LineWidth',1.3,'Color',[64/255 123/255 208/255]);
set(plot1(2),'DisplayName','Eq. (61)','MarkerSize',2,'LineWidth',2.3,...
    'LineStyle','-.',...
    'Color',[163/255 42/255 49/255]);

xlim([0,100])
xlabel('$\Gamma t$','Interpreter','latex')
ylabel('$|C_e(t)|^{2}$','Interpreter','latex')

% 
 %xlim(axes1,[0 50]);
% 
% ylim([0 1.05]);

% 
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
pbaspect([4 2.5 1])
% 
lgd = legend;
lgd.FontSize = 11;
lgd.Interpreter = 'latex';
legend('boxoff');
title('(d)','position',[6,0.36],'FontSize',14);

hold on
