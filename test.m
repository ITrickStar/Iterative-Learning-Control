%{
addpath(genpath('D:\packages\sedumi'));
addpath(genpath('D:\packages\YALMIP'));
%}
clear all; clc; close all;

Ts = 0.01;
Time = 1;
T = Time/Ts;
N = 100;
d = 3;

% график желаемой траектории
spx = 0:20:200;
spy = [0 0.01 0.042 0.1 0.14 0.15 0.14 0.1 0.042 0.01 0];
sp = spline(spx,spy);
xx = linspace (0,200,T+1);
yref = ppval(sp,xx)';

% параметры передаточной функции
W = tf([23.7356, 23.7356*661,2],[1,426.7,1.744*10^5,0]);
Y = ss(W);
sysd = c2d(Y,Ts);
[A, B, C] =  ssdata(sysd);

A_=[ A, zeros(3,3); -C*A, eye(1,3)];
B_ = [B; -C*B];
C_ = [C, zeros(1,3); zeros(3,3), eye(3)];

Q = diag([1 1 1 50]);
R = 10^-3;
K1 = -962.2;
K2 = 185.27;

%{
% матричное неравенство
Cond11 = X;
Cond12 = (X*A_+B_*Y*C_)';
Cond13 = X;
Cond14 = (Y*C_)';
Cond21 = A_*X+B_*Y*C_;
Cond22 = X;
Cond31 = X;
Cond33 = inv(Q);
Cond41 = YC_;
Cond44 = inv(R);
Cond = [Cond11, Cond12, Cond13, Cond14; ...
        Cond21, Cond22, 0,      0; ...
        Cond31, 0,      Cond33, 0 ...
        Cond41, 0,      0,      Cond44];
options = sdpsettings('solver', 'sedumi', 'verbose');
solved = solvesdp(Cond >= 0, options);
%}
%{
% начальные значения
x0 = 0;
x = zeros(1,N); xold = x;
E = zeros(1,N); Eold = E;

for t=1:T
    u(t,k)=u(t,k-1)+K(1)*(xold{t,k}-old{t,k-1})+K(2)*(yref(t+1,1)-C*x{t+1,k-1});
    x{t+1,k}=A*x{t,k}+B*u(t,k);
    y(t,k)=C*x{t,k};

    xold{1,t}=A*xold{t,k}+B*u(t,k)-C*xold{t,k};
    e(1,t-1)= yref(1,t)-y(1,t-1);
    E(1,t-1) = (E(1,t-1) + (e(t,k-1)^2));
end
    E(1,k-1) = sqrt(E(1,k-1)/N); % среднеквадратическая ошибка обучения
%}

E = zeros(1,N);
x0 = [0;0;0];
xhat0 = [0;0;0];

for k=1:N+1
    x{1,k}=x0;
    xhat{1,k}=xhat0;
    y(1,k)=C*x0;
end

for t=1:T+1

    u(t,1)=0;
    x{t+1,1}=A*x{t,1};
    xhat{t+1,1}=A*xhat{t,1};
    y(t,1)=C*x{t,1};
end

for k=2:N+1
    for t=1:T
        T_1=xhat{t,k-1};
        T_2=xhat{t,k};

        u(t,k)=u(t,k-1)+K1*(xhat{t,k}-xhat{t,k-1})+K2*(yref(t+1,1)-C*x{t+1,k-1});
        x{t+1,k}=A*x{t,k}+B*u(t,k);
        y(t,k)=C*x{t,k};
        xhat{t+1,k}=A*xhat{t,k}+B*u(t,k)+y(t,k)-C*xhat{t,k};
        e(t,k-1)= yref(t,1)-y(t,k-1);
        E(1,k-1) = (E(1,k-1) + (e(t,k-1)^2));          
    end
    E(1,k-1) = sqrt(E(1,k-1)/N);
end


% построение графика
[XX, YY] = meshgrid(0:T-1,0:N-1);
figure(1);
plot(YY, E', 'b','LineWidth', 2);
title("TEST", 'FontSize',14);
xlabel('k','FontSize',14);
ylabel('E','FontSize',14);
grid on;
