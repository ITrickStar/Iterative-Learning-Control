clear all; clc; close all;
Ts = 0.01;
Time = 2;
T = Time/Ts;
N = length(T);
d = 3;

x = 0:20:200;
y = [0 0.01 0.042 0.1 0.14 0.15 0.14 0.1 0.042 0.01 0];
cs = spline(x,y);
xx = linspace (0,T,T);
p = ppval(cs,xx);

plot(p);
xlabel('x','FontSize',15);
ylabel('y','FontSize',15);
grid on