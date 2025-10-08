% clear all
clc;
klist=0:8;
tau=1;  
omegae=0.5818*pi/tau;omegax=2*pi/tau;N=3; omegal=1.2*pi; Omega_3ls=-0.6*pi;Gamma=0.3*pi/tau; 
for j=1:length(klist)
    k=klist(j);
y1(j)=(omegax*tau-2*k*pi/N-tau*Omega_3ls^2/(omegae+omegal-2*k*pi/(N*tau)))/(-0.5*N*Gamma*tau);
end
figure(2)
subplot(2,2,1)
ezplot(@(k)cot(k*pi/N),[0 8])
hold on
plot(klist,y1,'o')
xlabel('$j$','interpreter','latex','fontsize',14);
hold on
%%
klist=10:40;
tau=1;  
omegae=0.2*pi/tau; omegax=7.0366*pi;  N=6; omegal=2.5*pi; Omega_3ls=2.8*pi;Gamma=0.1825*pi/tau;
for j=1:length(klist)
    k=klist(j);
y2(j)=(omegax*tau-2*k*pi/N-tau*Omega_3ls^2/(omegae+omegal-2*k*pi/(N*tau)))/(-0.5*N*Gamma*tau);
end
figure(2)
subplot(2,2,3)
ezplot(@(k)cot(k*pi/N),[15 40])
hold on
plot(klist,y2,'o')
hold on
xlabel('$j$','interpreter','latex','fontsize',14);