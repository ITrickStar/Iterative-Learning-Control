clear all; clc; close all;

% ������� ���������
Ts = 0.01;  % ���
Time = 0.1;  % �����
T = Time/Ts;    % ���������� ��������
N = 100;    % ���������� ����������
nx = 3;
d = 3;  % �������������

% ������ �������� ���������� �������� ������-�������������
spx = 0:20:200;
spy = [0 0.1 0.42 1 14 15 4 1 0.42 0.1 0];
sp = spline(spx,spy);
xx = linspace(0,T,T);
yref = ppval(sp,xx)';

% ��������� ������������ �������
W = tf([23.7356, 23.7356*661,2],[1,426.7,1.744*10^5,0]);
Y = ss(W);
sysd = c2d(Y,Ts);
[A, B, C] =  ssdata(sysd);

% ������� �������
Q = diag([1 1 1 50]);
R = 10^-3;
K = [-962.2, 185.27];

% �������� �� ������������
A_= [A, zeros(nx,1); -C*A, eye(1)];
B_ = [B; -C*B];
C_ = [C, 0; zeros(nx,nx), eye(3,1)];

X1 = sdpvar(nx,nx);
X2 = sdpvar(1);
X = [X1 zeros(nx,1); zeros(1,nx) X2];

% ��������� �����������
C11 = X;
C12 = (A_*X)';
C13 = X;
C21 = A_*X
C22 = X;
C31 = X;
C33 = inv(Q);

Cond = [C11,        C12,           C13;          ...
        C21,        C22,           zeros(4,4);  ...
        C31,        zeros(4,4),    C33;]

options = sdpsettings('solver', 'sedumi', 'verbose', 0);
solved = solvesdp(Cond >= 0, [], options);


% ��� ��������
E = zeros(1, N); % ������ ��������
x0 = zeros(nx,1);

% ������� ��������� ��������
for k=1:N+1
    x{1,k}=x0;
    y(1,k)=C*x0;
end

% �������� �� ������ ��������
for t=1:T
    u(t,1)=0;
    x{t+1,1}=A*x{t,1};
    y(t,1)=C*x{t,1};
end

% ��� ��������
% ��� ������ d-�����
for k=2:N+1
    for t=1:d
        u(t,k)=u(t,k-1)+K(1)*(C*x{t,k}-C*x{t,k-1})+K(2)*(yref(t+1,1)-C*x{t+1,k-1}); % ������ ����������
        x{t+1,k}=A*x{t,k}+B*u(t,k); % ��������� �������� �������
        y(t,k)=C*x{t,k};    % ��������� ��������� �������
        e(t,k-1)= yref(t,1)-y(t,k-1);   % ������ ��������������� �������
        E(1,k-1) = E(1,k-1) + (e(t,k-1)^2);     
    end
    E(1,k-1) = sqrt(E(1,k-1)/N); % ������ �������������������� ������
end

for k=2:N+1
    for t=d:T-1
        u(t,k)=u(t,k-1)+K(1)*(C*x{t-d,k}-C*x{t-d,k-1})+K(2)*(yref(t+1,1)-C*x{t-d+1,k-1}); % ������ ����������
        x{t+1,k}=A*x{t,k}+B*u(t,k); % ��������� �������� �������
        y(t,k)=C*x{t,k};    % ��������� ��������� �������
        e(t,k-1)= yref(t,1)-y(t,k-1);   % ������ ��������������� �������
        E(1,k-1) = E(1,k-1) + (e(t,k-1)^2);     
    end
    E(1,k-1) = sqrt(E(1,k-1)/N); % ������ �������������������� ������
end

% ���������� �������
[XX, YY] = meshgrid(0:T-1,0:N-1);
figure(1);
%plot(YY, y', 'b','LineWidth', 2)
plot(YY, E', 'b','LineWidth', 2);
title("ILC", 'FontSize',14);
xlabel('k','FontSize',14);
ylabel('E','FontSize',14);
grid on;
