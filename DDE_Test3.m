clc
clear all
IC = zeros(2*nc.n,1);
tspan = 1:0.001:10;                %Simulation time

    %sol    = dde23(    ddefile,    lags,        history,           tspan);
    sol      = dde23(       @dde,    [10],      @history,          tspan);

%% Plots of results  
ax1 = subplot(2,1,1);
    plot(sol.x,sol.y);
    title('Motion of each block of the system'); grid on; xlabel('Time t');    ylabel('x(t)');   
ax2 = subplot(2,1,2);
        xplot = 1:0.001:10;
        yplot = deval(sol(1:size(nc.n)),xplot,2*nc.n);
    plot(xplot,yplot);   grid on; xlabel('Time t');    ylabel('x(t)');
    title('Motion of the bit');

function s = history(t)
    s=zeros(2*nc.n,1);
    for i=1:2*nc.n
        s(i,1)=0;
    end
end

function xdot=dde(t,x,Z)                    % x - state of system (position of each block)
%% INPUT for DRILLSTRING    
% Stiffness matrix
% put last right bottom value as sum of spring (K) and "Cutting constant" (e)
K=full(gallery('tridiag',nc.n,-1,2,-1));    %spring
K(end,end)=1;                               %put last right bottom value to 1
K=nc.k*K;                                   %Stiffness matrix
K(end,end)=K(end,end)+nc.e;
%Eigenvalue decomposition
[V,D] = eig(inv(nc.M)*K);    
%% MVv0t                                 
v0t=zeros(nc.n,1);
v0t(1)=nc.v0*t;                                 % external force f1=f(t)*M^-1 - 9.81
VMv0t=inv(V)*inv(nc.M)*v0t;
%% VMVE
xlag=zeros(nc.n,1);
xlag(end,1)=Z(end,1);
VMVE=nc.e*inv(V)*inv(nc.M)*V*xlag;
%% VgB
gB=9.81*nc.B;
gBV=gB*inv(V);
gBV=diag(gBV);
%% D
D=D*eye(nc.n);
%% The differential equations:    
    x1      = zeros(nc.n,1);
    x2      = zeros(nc.n,1);
    x1dot   = zeros(nc.n,1);
    x2dot   = zeros(nc.n,1);
    xdot    = zeros(2*nc.n,1);
    x1    = x(1:nc.n);
    x2    = x(nc.n+1:2*nc.n);
    x1dot = (x2);
    x2dot = (VMv0t + VMVE - gBV - D*x1 -  nc.c*x2);
    xdot  = [x1dot; x2dot];
end