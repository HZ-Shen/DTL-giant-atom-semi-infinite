%
clear all
clc;
N=6; tau=1; 
k1=4*(N+1)-1;k2=4*(N+1)+2; 
tstep = 14000;                                                             
deltat = 1/tstep;                                                         
t = 80;                                                                                                                           
tlist = 0:deltat:t;   
Cx = zeros(size(tlist)); Ce = zeros(size(tlist)); 

%%%%%
x=1:1:N;
v=(x(2)-x(1))/tau;
X=0:0.01:N+1;
%%%%%

% Cx(1) = 1; Ce(1) = 0;
% omegax=2*pi/tau;k=1;                                                          
% Gamma=0.05*pi/tau;
% omegae=0;omegal=0; Omega_3ls=0;  

Cx(1) = sqrt(0.8); Ce(1) = sqrt(0.2); 
omegal=2.5*pi;Omega_3ls=2.8*pi;omegae=0.2*pi; 
omegax=Omega_3ls^2/(omegae+omegal-2*k1*pi/((N+1)*tau))-cot(k1*pi/(N+1))*(Omega_3ls^2/(omegae+omegal-2*k1*pi/((N+1)*tau))-Omega_3ls^2/(omegae+omegal-2*k2*pi/((N+1)*tau))+2*k1*pi/((N+1)*tau)-2*k2*pi/((N+1)*tau))/(cot(k1*pi/(N+1))-cot(k2*pi/(N+1)))+2*k1*pi/((N+1)*tau);
Gamma=2/(N+1)*(Omega_3ls^2/(omegae+omegal-2*k1*pi/((N+1)*tau))-Omega_3ls^2/(omegae+omegal-2*k2*pi/((N+1)*tau))+2*k1*pi/((N+1)*tau)-2*k2*pi/((N+1)*tau))/(cot(k1*pi/(N+1))-cot(k2*pi/(N+1)));
abs((cot(k2*pi/(N+1))*(N+1)*Gamma/2-Omega_3ls^2/(omegae+omegal-2*k2*pi/((N+1)*tau))-omegal/2)/omegax)                                                     
%%%%%
sj=Gamma*t;
sj1=fix(sj/10)*10;
[X1,Y1]=meshgrid(X,tlist);
p=zeros(size(X1));
R=1;
r=R+1i*sqrt(R*(1-R));
%%%%%

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
            index2 = j - (m+n)*tau*tstep;                                  %t-(m+n)*tau
            HH1=exp(1i*omegal/2*abs(m-n));HH2=exp(1i*omegal/2*(m+n));
            if tlist(j)-abs(m-n)*tau>0&&tlist(j)-(m+n)*tau<0               %|m-n|*tau<t<(m+n)*tau
                DeltaCx = DeltaCx - (Gamma)/2*HH1^tau*Cx(index1);
                
            elseif tlist(j)-(m+n)*tau>0                                    %t>(m+n)*tau
                DeltaCx = DeltaCx - (Gamma)/2*HH1^tau*Cx(index1)+ r*Gamma/2*HH2^tau*Cx(index2);
            
            else
                DeltaCx = DeltaCx;
            end
        end     
    end
   DeltaCx = DeltaCx + N*(Gamma)/2*Cx(j);  
   Cx(j+1)=Cx(j)+DeltaCx*deltat;                                           %b(t+deltat)=b(t)+Deltab*deltat
   Ce(j+1)=Ce(j)+DeltaCe*deltat;
end                                                                        %

for i=1:length(tlist)                                                      %
    t=tlist(i);                                                            %
for k=1:length(X)                                                          %
    sum1=0;                                                                %
    sum2=0;                                                                %
for m=1:N                                                                  %
    xm=x(m);                                                               %
    index1=t*tstep-abs(X(k)-xm)*tstep/v;                                   %
    index2=t*tstep-(X(k)+xm)*tstep/v; 
    HH1=exp(1i*omegal/2*abs(X(k)-xm));HH2=exp(1i*omegal/2*(X(k)+xm));
    if index1>0                                                           
        a1=fix(index1)+1;                                                  %                                             
        sum1=sum1+Cx(a1)*HH1;                                              %
    end
    if index2>0
        a2=fix(index2)+1;                                                  %                                            
        sum2=sum2+Cx(a2)*HH2;                                              %
    end
end
p(i,k)=Gamma/(2*v)*abs(-1i*sum1+1i*sum2)^2;                                %
end
end
p2=v*tau*p(1:70:end,:);
tlist1=tlist(1:70:end);
p2=p2(1:fix(sj1/sj*length(tlist1)),:);
Y1=Y1(1:fix(sj1/sj*length(tlist1)),:);
figure(1)
subplot('position',[0.12 0.1 0.366471449487555 0.378846153846154])
surf(p2)
shading interp
colormap(othercolor('BuDOr_12'))
xlim([0+0.5,size(Y1,2)])
ylim([0+0.5,size(Y1,1)+0.5])
zlim([0,0.4])
xtick=linspace(1,length(X),N-(-1)+1);
ytick=linspace(1,fix(sj1/sj*length(tlist1)),3);
caxis([0,0.4])
%axis xy
xticks(xtick)
yticks(ytick)
xticklabels(0:1:N+1)
yticklabels(linspace(0,sj1,3))
box on
grid on

xlabel('$x/x_{0}$','FontSize',12,'Interpreter','latex')
ylabel('$t/\tau_0$','FontSize',12,'Interpreter','latex')
zlabel('$x_{0}P(x,t)$','FontSize',12,'Interpreter','latex')
view([135,25])
hold on
title('(c)','position',[575.7752965214568,-93.04113690951272,0.37581614963482],'FontSize',14);

hold on
