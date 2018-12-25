clc
clear all
tspan = 1:1:100;
n=3;
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
    xplot = 1:.001:100;
    yplot = deval(sol,xplot,2*n);
    plot(xplot,yplot);
    xlabel('time t');
    ylabel('x(t)');

function xdot=dde(t,x,Z)                % x - state of system (position of each block)
%% INPUT for BIT
z_b=5;                                  %number of blades
a=0.4445/2;                             %radius of section#3 in [m]
                                        %Orientation of cutting face =1 at 1st iteration
se=1;                                   %MSE Intrinsic Energy SHOULD BE 20 NEED TO CHECK???????????
e=z_b*a*se;                             %Cutting Energy
%% INPUT for DRILLSTRING
v0=1;                                   %constant speed of drilling
c=0.2;                                  %damping coefficient
n=3;                                    %number of elements for section #3 - length 1350m       
m=1;                                    %mass of each block
M=m*eye(n)                              %Matrix of drillstring mass
k=1;
K=full(gallery('tridiag',n,-1,2,-1));   %spring
K=k*K;
K(end,end)=k+e;                         %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
                                        %spring of each block
                                        %Stiffness matrix
[V,D] = eig(K,M);                       %Eigenvalue decomposition
%% MVF1                                 
mv=inv(M)*inv(V);
mvzeros=zeros(n,n);
mvzeros(1,1)=mv(1,1);
F1=zeros(n,1); 
F1(1)=v0*t;
MVF1=mvzeros*F1
%% C
C=c*eye(n)
%% D
D=D*eye(n)
%% MExlag                           Energy for cutting - Wc only for last element of drillstring (n)
E=zeros(n,n);
E(end,end)=e;                           %last element of vector E is cutting constant e                                       
ME=inv(M)*E;
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
        x2dot = MVF1 - C*x2 -D*x1 + ME*xlag;
        xdot = [x1dot; x2dot];
end

function s = history(t)
n=3;
    s=zeros(2*n,1);
    for i=1:2*n
        s(i,1)=1;
    end
end