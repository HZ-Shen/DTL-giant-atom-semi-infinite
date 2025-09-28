%1个解的情况:fig3a (2kpi)
clear all
clc;
N=3; tau=1; 
tstep = 1000;                                                             
deltat = 1/tstep;                                                         
t = 70;                                                                                                                           
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
omegax=2*pi/tau;k=4;
omegal=1.2*pi;Omega_3ls=-0.6*pi;Gamma=0.3*pi/tau; 
omegae=2*Omega_3ls^2/(Gamma*N*cot(k*pi/N)-4*k*pi/(N*tau)+2*omegax)-omegal+2*k*pi/(N*tau);%omegae has expression                                                    

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
            HH1=exp(1i*omegal/2*abs(m-n));
            if tlist(j)-abs(m-n)*tau>0                                     %|m-n|*tau<t<(m+n)*tau
                DeltaCx = DeltaCx - (Gamma)/2*HH1^tau*Cx(index1);                            
            else
                DeltaCx = DeltaCx;
            end
        end     
    end
   DeltaCx = DeltaCx + N*(Gamma)/2*Cx(j);  
   Cx(j+1)=Cx(j)+DeltaCx*deltat;                                           %b(t+deltat)=b(t)+Deltab*deltat
   Ce(j+1)=Ce(j)+DeltaCe*deltat;
end                                                                        %模拟β的数组

for i=1:length(tlist)                                                      %纵坐标时间变化的for循环
    t=tlist(i);                                                            %某一时刻
for k=1:length(X)                                                          %横坐标位置变化的for循环
    sum1=0;                                                                %|X-Xm|求和项总和
for m=1:N                                                                  %求和项中每一项的循环
    xm=x(m);                                                               %第m个耦合点
    index1=t*tstep-abs(X(k)-xm)*tstep/v;                                   %求和项中Beta数组索引
    HH1=exp(1i*omegal/2*abs(X(k)-xm)/v);
    if index1>0                                                           
        a1=fix(index1)+1;                                                  %索引取整数+1                                                 
        sum1=sum1+Cx(a1)*HH1;                                              %利用循环求和
    end
end
p(i,k)=Gamma/(2*v)*abs(-1i*sum1)^2;                                        %按照横坐标和纵坐标for循环分别向p矩阵填入数据
end
end
%变换坐标轴
%plot(tlist,abs(Cx).^2)
p2=v*tau*p(1:70:end,:);
tArray1=tlist(1:70:end);
p2=p2(1:fix(sj1/sj*length(tArray1)),:);
Y1=Y1(1:fix(sj1/sj*length(tArray1)),:);
p_2d=p2(end,:);
figure1=figure(1);
axes2 = axes('Position',[0.60 0.16 0.366471449487555 0.378846153846154]);
area(X/v/tau,p_2d,...
    'FaceColor',[205/255 205/255 193/255],...
    'EdgeColor',[139/255 139/255 131/255],...
    'LineWidth',2,'LineStyle','-.');
ylim(axes2,[0 0.4])
set(axes2,'FontName','Times New Roman','FontSize',12,'Layer','top',...
    'LineWidth',1,'XTick',[0 1 2 3 4],'YTick',[0 0.1 0.2 0.3 0.4])
%pbaspect([4 2.5 1])
xlabel('$x/x_0$','Interpreter','latex')
ylabel('$x_0 P(x,t\rightarrow \infty)$','Interpreter','latex')
title('(b)','position',[2.5,0.3],'FontSize',14);