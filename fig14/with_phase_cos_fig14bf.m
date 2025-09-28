clc;
clear;
close all;
tic;

% 程序中所提到的方程Eq .(10)和Eq .(11)为论文---AX12404孙思杰的论文250103.pdf---里的。

% parameters setting
tau_0 = 1;
N = 3;
j = 4;
G = -0.6*pi/tau_0;
Gamma = 0.3*pi/tau_0;
omega_x = 2*pi/tau_0;
omega_l = 6.46603;
omega_e = 0.5818*pi;
mu_1 = omega_x-0.5*omega_l;
mu_2 = omega_e+0.5*omega_l;
alpha = 0;
t_end = 300;


% Define time span for simulation in units of t (will scale later by Gamma)
tspan = [0, t_end/tau_0]; 

% Define initial condition (history function) 
history = @(t) [sqrt(0.8); sqrt(0.2)];

% Solve the DDE using ddesd solver
options = ddeset(RelTol=1e-3);  % adjust the tolerance when needed!
sol = ddesd(ddefunf(tau_0, Gamma, mu_1, mu_2, G, omega_l, N, alpha), ddelagsf(tau_0, N), history, tspan, options);

tf = t_end/tau_0;
t = linspace(0, tf, 1000);
y = deval(sol, t);

% Calculate the atomic excitation
Final_x = y(1,:).*conj(y(1,:));   %Cxt的纯数值解
Final_e = y(2,:).*conj(y(2,:));   %Cet的纯数值解



% 渐进稳态解Cxt
eps=zeros(size(t)); 
omega_j = 6.19176; %对应于论文中的omega_j
s_j = -1i*omega_j;
cx0 = sqrt(0.8);
ce0 = sqrt(0.2);
constant = cxtfun(tau_0,N,G,Gamma,omega_l,mu_1,mu_2,alpha,cx0,ce0,s_j);
for i=1:length(t)
    eps(i)=constant;
end  
Final_xeq = eps;

% 渐进稳态解Cet
eps=zeros(size(t)); 
omega_j = 6.19176; %对应于论文中的omega_j
s_j = -1i*omega_j;
cx0 = sqrt(0.8);
ce0 = sqrt(0.2);
constant = cetfun(tau_0,N,G,Gamma,omega_l,mu_1,mu_2,alpha,cx0,ce0,s_j);
for i=1:length(t)
    eps(i)=constant;
end  
Final_eeq = eps;

% Plot the atomic excitation over t./tau_0
figure;
indxlist = t./tau_0;
plot(indxlist, Final_x, LineStyle="-", Color=[0.07,0.62,1.00], LineWidth=0.7)
hold on
plot(indxlist, Final_xeq, LineStyle="--", Color=[1.00,0.41,0.16], LineWidth=1)
xlim([0 300])
ylim([0 0.8])
xlabel('$t/\tau_0$', Interpreter='latex')
ylabel('$|C_{x}(t)|^2$', Interpreter='latex')
title('Fig. 3a (with phase)', Interpreter='latex');
legend('$Eq. (8)$','$C_{x}(t)$ steady-state solution', Interpreter='latex')
grid on;

figure;
indxlist = t./tau_0;
plot(indxlist, Final_e, LineStyle="-", Color=[0.07,0.62,1.00], LineWidth=0.7)
hold on
plot(indxlist, Final_eeq, LineStyle="--", Color=[1.00,0.41,0.16], LineWidth=1)
xlim([0 300])
ylim([0 0.8])
xlabel('$t/\tau_0$', Interpreter='latex')
ylabel('$|C_{e}(t)|^2$', Interpreter='latex')
title('Fig. 3a (with phase)', Interpreter='latex');
legend('$Eq. (9)$','$C_{e}(t)$ steady-state solution', Interpreter='latex')
grid on;
hold on

toc





%%%%%%%%%%%%%%%%%%%%%setting the deday differential equation%%%%%%%%%%%%%
% 当前延迟微分方程对应于 ---AX12404孙思杰的论文250103.pdf--- 中的Eq. (10)和Eq. (11)，注意加入了相位
% 如果延迟微分方程变了，只需要理解并修改这部分内容即可继续使用此程序 !!!
function dxdtf = ddefunf(tau_0, Gamma, mu_1, mu_2, G, omega_l, N, alpha)

dxdtf = @ddefun;

function dxdt = ddefun(t,epsilon,Z)
%%%%%%%%
s1 = 0;
s2 = 0;
s3 = 0;
s4 = 0;

for m = 1:N
    t1 = -(Gamma/2).*epsilon(1);
    s1 = s1+t1;
end

for m = 2:N
    for n = 1:m-1
        count = (m-1)*(m-2)/2;       % count 的作用是跳过当前m之前的所有已赋值的数组元素
        x = n+count;                 % x代表了与m-n相关的延迟项在Z这个数组中的位置，这个位置(Z中的第二个指标)与相应的与m-n相关的延迟项所对应的延迟是什么有关，(N^2-N)/2是延迟项m-n(或m+n)所需要的数组大小
        y = x+((N^2-N)/2);           % y代表了与m+n相关的延迟项在Z这个数组中的位置，这个位置(Z中的第二个指标)与相应的与m+n相关的延迟项所对应的延迟是什么有关
%         t2 = -(Gamma).*cos(alpha.*n-alpha.*m).*Z(1,x).*heaviside(t-abs(m-n).*tau_0).*exp(1i.*0.5.*omega_l.*abs(m-n).*tau_0);
%         t4 = (Gamma).*cos(alpha.*n-alpha.*m).*Z(1,y).*heaviside(t-(m+n).*tau_0).*exp(1i.*0.5.*omega_l.*(m+n).*tau_0);
        t2 = 0;
        t4 = 0;
        s2 = s2+t2;
        s4 = s4+t4;
    end
end

for m = 1:N 
    y1 = (N^2-N)+m;   % y1代表了与2*m相关的延迟项在Z这个数组中的位置，这个位置(Z中的第二个指标)与相应的与2*m相关的延迟项所对应的延迟是什么有关，(N^2-N)是延迟项m-n加延迟项m+n一共所需要的数组大小
    t3 = (Gamma/2).*Z(1,y1).*heaviside(t-(2*m).*tau_0).*exp(1i.*0.5.*omega_l.*(2*m).*tau_0);
    s3 = s3+t3;
end

dxdt = [-1i*mu_1*epsilon(1)-1i*G*epsilon(2)+s1+s2+s3+s4;
        -1i*mu_2*epsilon(2)-1i*G*epsilon(1)];
end

end




%%%%%%%%%%%%% non constant delay%%%%%%%%%%%%%%%%%%%
% 如果延迟微分方程变了，只需要理解并修改这部分内容即可继续使用此程序 !!!
function delayf = ddelagsf(tau_0,N)

delayf = @ddelags;

function delay = ddelags(t,epsilon)
delay = zeros([1 N^2]);
for m = 2:N
    for n = 1:m-1
        count = (m-1)*(m-2)/2;
        x = n+count;   % 这里对应在数组Z位置x处的延迟项
        y = x+((N^2-N)/2);           % 这里对应在数组Z位置y处的延迟项
        delay(1,x) = t-abs(m-n)*tau_0;   % 赋值对应于数组Z在位置x处的延迟项的延迟
        delay(1,y) = t-(m+n)*tau_0;      % 赋值对应于数组Z在位置y处的延迟项的延迟
    end
end

for m = (N^2-N)+1:(N^2-N)+N
        y1 = m;                           % 这里对应在数组Z位置y1处的延迟项
        m1 = m-(N^2-N);
        delay(1,y1) = t-(2*m1)*tau_0;      % 赋值对应于数组Z在位置y1处的延迟项的延迟
end

delay = delay';
end

end

%%%%%%%%%%%Cxt渐进稳态解%%%%%%%%%%%%%
function result = cxtfun(tau0,N1,G,gamma,omegal,mu1,mu2,alpha,cx0,ce0,s)
    x1 = (sqrt(-1)*(-1)).*ce0.*G+cx0.*(sqrt(-1).*mu2+s);                      %对应于Eq. 13的分子
    us1 = (1/2).*exp(1).^((-2).*s.*tau0+(-2).*N1.*s.*tau0).*((-1)+exp(1).^( ...
  2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0)).^(-1).*((-1).*exp(1) ...
  .^(sqrt(-1).*omegal.*tau0+sqrt(-1).*N1.*omegal.*tau0)+exp(1).^(2.* ...
  ((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0+2.*s.*tau0+2.*N1.*s.* ...
  tau0)).*gamma+sqrt(-1).*mu1+(1/2).*gamma.*N1+s+(sqrt(-1).*mu2+s).* ...
  (1+exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0+(-2).*s.* ...
  tau0+(-2).*N1.*s.*tau0).*((-1)+exp(1).^(2.*((sqrt(-1)*(1/2)).* ...
  omegal+(-1).*s).*tau0)).^(-2).*((-1).*exp(1).^(sqrt(-1).*omegal.* ...
  tau0+sqrt(-1).*N1.*omegal.*tau0)+exp(1).^(2.*((sqrt(-1)*(1/2)).* ...
  omegal+(-1).*s).*tau0+2.*s.*tau0+2.*N1.*s.*tau0)).*gamma.*tau0+ ...
  exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0).*((-1)+exp( ...
  1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0)).^(-1).*gamma.* ...
  N1.*tau0+(1/2).*exp(1).^((-2).*s.*tau0+(-2).*N1.*s.*tau0).*((-1)+ ...
  exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0)).^(-1).*(( ...
  -1).*exp(1).^(sqrt(-1).*omegal.*tau0+sqrt(-1).*N1.*omegal.*tau0)+ ...
  exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0+2.*s.*tau0+ ...
  2.*N1.*s.*tau0)).*gamma.*((-2).*tau0+(-2).*N1.*tau0));
    result = abs(x1./us1)^2;
end

%%%%%%%%%%%Cet渐进稳态解%%%%%%%%%%%%%
function result = cetfun(tau0,N1,G,gamma,omegal,mu1,mu2,alpha,cx0,ce0,s)
    x1 = (sqrt(-1)*(-1)).*cx0.*G+ce0.*((1/2).*exp(1).^((-2).*s.*tau0+(-2).* ...
  N1.*s.*tau0).*((-1)+exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s) ...
  .*tau0)).^(-1).*((-1).*exp(1).^(sqrt(-1).*omegal.*tau0+sqrt(-1).* ...
  N1.*omegal.*tau0)+exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).* ...
  tau0+2.*s.*tau0+2.*N1.*s.*tau0)).*gamma+sqrt(-1).*mu1+(1/2).* ...
  gamma.*N1+s);
    us1 = (1/2).*exp(1).^((-2).*s.*tau0+(-2).*N1.*s.*tau0).*((-1)+exp(1).^( ...
  2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0)).^(-1).*((-1).*exp(1) ...
  .^(sqrt(-1).*omegal.*tau0+sqrt(-1).*N1.*omegal.*tau0)+exp(1).^(2.* ...
  ((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0+2.*s.*tau0+2.*N1.*s.* ...
  tau0)).*gamma+sqrt(-1).*mu1+(1/2).*gamma.*N1+s+(sqrt(-1).*mu2+s).* ...
  (1+exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0+(-2).*s.* ...
  tau0+(-2).*N1.*s.*tau0).*((-1)+exp(1).^(2.*((sqrt(-1)*(1/2)).* ...
  omegal+(-1).*s).*tau0)).^(-2).*((-1).*exp(1).^(sqrt(-1).*omegal.* ...
  tau0+sqrt(-1).*N1.*omegal.*tau0)+exp(1).^(2.*((sqrt(-1)*(1/2)).* ...
  omegal+(-1).*s).*tau0+2.*s.*tau0+2.*N1.*s.*tau0)).*gamma.*tau0+ ...
  exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0).*((-1)+exp( ...
  1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0)).^(-1).*gamma.* ...
  N1.*tau0+(1/2).*exp(1).^((-2).*s.*tau0+(-2).*N1.*s.*tau0).*((-1)+ ...
  exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0)).^(-1).*(( ...
  -1).*exp(1).^(sqrt(-1).*omegal.*tau0+sqrt(-1).*N1.*omegal.*tau0)+ ...
  exp(1).^(2.*((sqrt(-1)*(1/2)).*omegal+(-1).*s).*tau0+2.*s.*tau0+ ...
  2.*N1.*s.*tau0)).*gamma.*((-2).*tau0+(-2).*N1.*tau0));
    result = abs(x1./us1)^2;
end
