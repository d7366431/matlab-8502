clc
clear all
IC = zeros(2*nc.n,1);
tspan = 1:0.01:4;                %Simulation time

    %sol    = dde23(    ddefile,    lags,        history,           tspan);
    sol      = dde23(       @dde,    [4],      @history,          tspan);
    [t,x_nd] = ode45(    @dde_nd,     tspan,      IC);

%% Plots of results  
ax1 = subplot(3,1,1);
    plot(sol.x,sol.y);
    title('sol.x, sol.y'); grid on; xlabel('Time t');    ylabel('x(t)');   
ax2 = subplot(3,1,2);
        xplot = 1:0.01:4;
        yplot = deval(sol(1:size(nc.n)),xplot,2*nc.n);
    title('deval');
    plot(xplot,yplot);   grid on; xlabel('Time t');    ylabel('x(t)');
ax2 = subplot(3,1,3);
    plot(t,x_nd(:,end)); grid on; xlabel('Time t no delay'); ylabel('x(t) no delay');

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
v0t=nc.v0*t;                                 % external force f1=f(t)*M^-1 - 9.81
MVv0t=zeros(nc.n,1);
MVv0t(1)=v0t;
MVv0t=inv(nc.M)*inv(V)*MVv0t;
%% VMVE without lag                          Energy for cutting - Wc only for last element of drillstring (n)
VMVE=nc.e*inv(V)*inv(nc.M)*V;
%% VgB
gB=9.81*nc.B;
VGB=inv(V)*gB;
VGB=diag(VGB);
%% D
D=D*eye(nc.n);
%% The differential equations:
    xlag    = Z(1:nc.n,1);    
    x1      = zeros(nc.n,1);
    x2      = zeros(nc.n,1);
    x1dot   = zeros(nc.n,1);
    x2dot   = zeros(nc.n,1);
    xdot    = zeros(2*nc.n,1);
    x1 = x(1:nc.n);
    x2 = x(nc.n+1:2*nc.n);
    x1dot = x2;
    x2dot = MVv0t + VMVE*xlag - D*x1 -  nc.c*x2 - VGB;
    %x2dot_nd        = MVv0t + VMVE*x1_nd_zeros %- D*x1_nd -  nc.C*x2_nd;% - VGB;       %Lag is zero %%Missing VGB - GETTING ERROR
    xdot            = [x1dot; x2dot];
end

function xdot_nd=dde_nd(t,x_nd)             % x - state of system (position of each block)
%% INPUT for DRILLSTRING
% Stiffness matrix
%put last right bottom value as sum of spring (K) and "Cutting constant" (e)
K=full(gallery('tridiag',nc.n,-1,2,-1));    %spring
K(end,end)=1;                               %put last right bottom value to 1
K=nc.k*K;                                   %Stiffness matrix
K(end,end)=K(end,end)+nc.e;
%Eigenvalue decomposition
[V,D] = eig(inv(nc.M)*K);    
%% MVv0t                                 
v0t=nc.v0*t;                                 % external force f1=f(t)*M^-1 - 9.81
MVv0t=zeros(nc.n,1);
MVv0t(1)=v0t;
MVv0t=inv(nc.M)*inv(V)*MVv0t;
%% VMVE without lag                          Energy for cutting - Wc only for last element of drillstring (n)
VMVE=nc.e*inv(V)*inv(nc.M)*V;
%% VgB
gB=9.81*nc.B;
VGB=inv(V)*gB;
VGB=diag(VGB);
%% D
D=D*eye(nc.n);
%VDV=V*D*inv(V);
%% The differential equations:
    x1_nd               = zeros(nc.n,1);
    x2_nd               = zeros(nc.n,1);
    x1_nd_zeros         = zeros(nc.n,1);
    x1dot_nd            = zeros(nc.n,1);
    x2dot_nd            = zeros(nc.n,1);
    x1dot_nd_zeros      = zeros(nc.n,1);
    xdot_nd             = zeros(2*nc.n,1);
    x1_nd               = x_nd(1:nc.n);
    x2_nd               = x_nd(nc.n+1:2*nc.n);
    x1_nd_zeros(end,1)  = x1_nd_zeros(nc.n,1);
    x1dot_nd            = x2_nd;
    x2dot_nd            = MVv0t + VMVE*x1_nd_zeros - D*x1_nd -  nc.c*x2_nd - VGB;
    %x2dot_nd           = MVv0t + VMVE*x1_nd_zeros %- D*x1_nd -  nc.C*x2_nd;% - VGB;       %Lag is zero %%Missing VGB - GETTING ERROR
    xdot_nd             = [x1dot_nd; x2dot_nd];
end

function s = history(t)
    s=zeros(2*nc.n,1);
    for i=1:2*nc.n
        s(i,1)=0;
    end
end