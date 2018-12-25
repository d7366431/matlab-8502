clc
clear all
tspan = 1:1:15;
n=30;
%%sol   = dde23(    ddefile,    lags,        history,           tspan);
  sol   = dde23(       @dde,    [1],        @history,           tspan);

  
ax1 = subplot(2,1,1);
plot(sol.x,sol.y);
grid on;
title('sol.x, sol.y')
    xlabel('time t');
    ylabel('x(t)');
 
ax2 = subplot(2,1,2);
grid on;
title('deval')
    xplot = 1:1:10;
    yplot = deval(sol,xplot,2*n);
    plot(xplot,yplot);
    xlabel('time t');
    ylabel('x(t)');

function xdot=dde(t,x,Z)                % x - state of system (position of each block)
n=30; %number of elements for simplicity put 3 for now       
v0=zeros(n,1);
v0(1)=1;
% c matrix
c=1;
C=zeros(n);
C(end,end)=c;
% D Matrix
d=1;
D=d*eye(n);
% M Matrix
m=1;                                   
M=m*eye(n);
%% The differential equations:
        xlag    = Z(1:n,1);
        x1      = zeros(n,1);
        x2      = zeros(n,1);
        x1dot   = zeros(n,1);
        x2dot   = zeros(n,1);
        xdot    = zeros(2*n,1);
        x1 = xdot(1:n);
        x2 = xdot(n+1:end);
        x1dot = x2;
        x2dot = v0*t - C*x2 -D*x1 - M*xlag; %v0*t-M*xlag-D*x1-C*x2;
        xdot = [x1dot; x2dot];
end

function s = history(t)
n=30;
    s=zeros(2*n,1);
    for i=1:2*n
        s(i,1)=1;
    end
end