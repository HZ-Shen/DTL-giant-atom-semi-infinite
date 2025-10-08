clear all
clc;
klist=0:8;
tau=1;  
omegae=0.5222*pi/tau;omegax=1.2*pi/tau;N=3; omegal=0.7*pi; Omega_3ls=0.5*pi;Gamma=0.3*pi/tau; 
for j=1:length(klist)
    k=klist(j);
y1(j)=(omegax*tau-2*k*pi/(N+1)-tau*Omega_3ls^2/(omegae+omegal-2*k*pi/((N+1)*tau)))/(-0.5*(N+1)*Gamma*tau);
end
figure(2)
subplot(2,2,2)
ezplot(@(k)cot(k*pi/(N+1)),[0 8])
hold on
plot(klist,y1,'o')
xlabel('$j$','interpreter','latex','fontsize',14);
hold on
%%
clear all
clc;
klist=10:40;
tau=1;  
omegae=0.2*pi/tau; omegax=6.9350*pi; N=6; omegal=2.5*pi; Omega_3ls=2.8*pi;Gamma=0.1079*pi/tau;
for j=1:length(klist)
    k=klist(j);
y1(j)=(omegax*tau-2*k*pi/(N+1)-tau*Omega_3ls^2/(omegae+omegal-2*k*pi/((N+1)*tau)))/(-0.5*(N+1)*Gamma*tau);
end
figure(2)
subplot(2,2,4)
ezplot(@(k)cot(k*pi/(N+1)),[15 40])
hold on
plot(klist,y1,'o')
xlabel('$j$','interpreter','latex','fontsize',14);
hold on

