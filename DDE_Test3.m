clc
clear all
IC = zeros(2*nc.n,1);
tspan = 1:0.01:40;                %Simulation time

%%sol   = dde23(    ddefile,    lags,        history,           tspan);
  sol   = dde23(       @dde,    [0.01],        @history,           tspan);

[t,x_nd] = ode45(@dde_nd,tspan,IC);

%% Plots of results  
ax1 = subplot(3,1,1);
        plot(sol.x,sol.y);
    grid on;
    title('sol.x, sol.y')
    xlabel('Time t');    ylabel('x(t)');
 
ax2 = subplot(3,1,2);
        xplot = 1:0.01:40;
        yplot = deval(sol(1:size(nc.n)),xplot,2*nc.n);
    grid on;
    title('deval')
    plot(xplot,yplot);
    xlabel('Time t');    ylabel('x(t)');

ax2 = subplot(3,1,3);
    plot(t,x_nd(:,end)); 
    xlabel('Time t no delay'); ylabel('x(t) no delay');   
    
function xdot=dde(t,x,Z)                   % x - state of system (position of each block)
%% INPUT for DRILLSTRING    
%Stiffness matrix
K=full(gallery('tridiag',nc.n,-1,2,-1));   %spring
K(end,end)=1;
K=nc.k*K;
K(end,end)=K(end,end)+nc.e;                %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
%Eigenvalue decomposition
[V,D] = eig(K,nc.M);    
%% MVF1      
f1=nc.v0*t;                                % external force f1=f(t)*M^-1 - 9.81
F1=zeros(nc.n,1);
F1(1)=f1;
MVF1=inv(nc.M)*inv(V)*F1;
%% D
D=D*eye(nc.n);
%% MExlag                               Energy for cutting - Wc only for last element of drillstring (n)
E=zeros(nc.n,nc.n);
E(end,end)=nc.e;                           %last element of vector E is cutting constant e                                       
ME=inv(nc.M)*E;
%% The differential equations:
        xlag    = Z(1:nc.n,1);
        x1      = zeros(nc.n,1);
        x2      = zeros(nc.n,1);
        x1dot   = zeros(nc.n,1);
        x2dot   = zeros(nc.n,1);
        xdot    = zeros(2*nc.n,1);
        x1 = xdot(1:nc.n);
        x2 = xdot(nc.n+1:end);
        x1dot = x2;
        x2dot = MVF1 - nc.C*x2 -D*x1 + ME*xlag;
        xdot = [x1dot; x2dot]
        
end

function xdot_nd=dde_nd(t,x_nd)                % x - state of system (position of each block)
%% INPUT for DRILLSTRING
%Stiffness matrix
              %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
K=full(gallery('tridiag',nc.n,-1,2,-1));   %spring
K=nc.k*K;
K(end,end)=nc.k+nc.e;
%Eigenvalue decomposition
[V,D] = eig(K,nc.M);    
%% MVF1                                 
f1=nc.v0*t;                                % external force f1=f(t)*M^-1 - 9.81
F1=zeros(nc.n,1);
F1(1)=f1;
MVF1=inv(nc.M)*inv(V)*F1;

%% D
D=D*eye(nc.n);
%% The differential equations:
    x1_nd      = zeros(nc.n,1);
    x2_nd      = zeros(nc.n,1);
    x1dot_nd   = zeros(nc.n,1);
    x2dot_nd   = zeros(nc.n,1);
    xdot_nd    = zeros(2*nc.n,1);
    x1_nd = x_nd(1:nc.n);
    x2_nd = x_nd(nc.n+1:end);
    x1dot_nd = x2_nd;
    x2dot_nd = MVF1-nc.C*x2_nd - D*x1_nd; %Lag is zero
    xdot_nd = [x1dot_nd; x2dot_nd];
end

function s = history(t)
    s=zeros(2*nc.n,1);
    for i=1:2*nc.n
        s(i,1)=1;
    end
end