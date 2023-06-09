clear all; clc; close all;

Ts = 0.01;
Time = 3; 
T = Time/Ts; 
N = 100; 
nx = 4;

C_r = 0.1;
C_l = 0.1;
B_eq = 0.004;

Ac = [0 0 1 0; 0 0 0 1; 0 623.7*C_r -7.20*B_eq 0; 0 -965.3*C_l 7.20*B_eq 0];
Bc = [0; 0; 479.8052; -479.8052];
Cc = [1 0 0 0];
Dc = zeros(1,1);
%?????????? ??????
sys = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sys,Ts);
A = sysd.a;
B = sysd.b;
C = sysd.c;
D = sysd.d;

%??????? ???????
Qn=1e-2; 
Rn=1e-6; 

Dn = Bc;

Qn = Dn * Qn *Dn';
Qn = (Qn+Qn')/2;
Rn = (Rn+Rn')/2;


M = [-Ac  Qn ; zeros(nx) Ac'];
phi = expm(M*Ts);
phi12 = phi(1:nx,nx+1:2*nx);
phi22 = phi(nx+1:2*nx,nx+1:2*nx);
Qd = phi22'*phi12;
Qd = (Qd+Qd')/2; 
Rd = Rn/Ts; 

Pstar =[2, 0, 0, 0;...
        0, 1, 0, 0;...
        0, 0, 3, 0;...
        0, 0, 0, 4];

for i=1:1000
  F = Pstar*C'/(C*Pstar*C'+Rd);
  Ptiled = Pstar-F*C*Pstar;
  Pstar = A*Ptiled*A'+Qd;
end


A11 = [A-F*C zeros(nx); F*C A];
A12 = zeros(8,1);
A21 = [-C*A -C*A];
A22 = ones(1,1);

A_ = [A11 A12; A21 A22];

B1 = [zeros(4,1); B];
B2 = -C*B;
B_ = [B1;B2];


X1 = sdpvar(2*nx,2*nx);
X2 = sdpvar(9-2*nx,9-2*nx);
X = [X1 zeros(2*nx,9-2*nx); zeros(9-2*nx,2*nx) X2];

sigma = 0.5;

q = [1e-1 2e-1 5e-1 2e-1 6e-1 5e-1 2e-1 3e-1 100];

Q =diag(q);
R = 1e-3;

S11 = (1-sigma)*X;
S12 = X*A_';
S13 = X;
S21 = A_*X;
S22 = X + B_*inv(R)*B_';
S23 = zeros(9,9);
S31 = X;
S32 = zeros(9,9);
S33 = inv(Q);

Og = [[S11, S12, S13; S21, S22, S23; S31, S32, S33]>=0, X>=1e-16*eye(9)] ;

options = sdpsettings('solver', 'sedumi', 'verbose', 0);
diagn = solvesdp(Og,[],options);
disp(diagn);

Res = 0;

if diagn.problem == 0
   Res = double(X);
end

P=inv(Res);

spx = 0:20:200;
spy = [0 0.01 0.042 0.1 0.14 0.15 0.14 0.1 0.042 0.01 0];
sp = spline(spx,spy);
xx = linspace (0,200,T+1);
yref = ppval(sp,xx)';


E = zeros(1,N);

teta3=eye(1);

K = (-1)*inv(B_'*P*B_+R)* B_'*P*A_;
K2=K(1, [5, 6, 7, 8]);
K3=K(1,9);


x0 = [0;0;0;0];
xhat0 = [0;0;0;0];

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

%?????????????
for k=2:N+1
    for t=1:T
%????????
        T_1=xhat{t,k-1};
        T_2=xhat{t,k};

        u(t,k)=u(t,k-1)+K2*(xhat{t,k}-xhat{t,k-1})+K3*teta3*(yref(t+1,1)-C*x{t+1,k-1});
%?????? ???????
        x{t+1,k}=A*x{t,k}+B*u(t,k);
        y(t,k)=C*x{t,k};
%????????? ??????
        xhat{t+1,k}=A*xhat{t,k}+B*u(t,k)+F*(y(t,k)-C*xhat{t,k});
%?????? ??????????????? ???????
        e(t,k-1)= yref(t,1)-y(t,k-1);
%??????????? ??????
        E(1,k-1) = (E(1,k-1) + (e(t,k-1)^2));          
    end
    E(1,k-1) = sqrt(E(1,k-1)/T);
end

%???????     
[XX, YY] = meshgrid(0:T-1,0:N-1);

figure(1);

plot(YY, E', 'b', 'LineWidth', 2);
xlabel('k','FontSize',15);
ylabel('E','FontSize',15);
grid on;

for k=1:N
   for t=1:T
       uu(t,k)=u(t,k);
       ee(t,k)=e(t,k);
   end
end

figure(2);

subplot(1,2,1);
mesh(XX, YY, uu' );
title('Change of control depending on the number of repetitions','FontSize', 15);
xlabel('t','FontSize', 15);
ylabel('k','FontSize', 15);
zlabel('u','FontSize', 15);
