clc
clear all
n=3;                               %number of elements for section #3 - length 1350m  
IC = zeros(2*n,1)
tspan = 1:0.001:40;                    %Simulation time
                                   %number of elements for section #3 - length 1350m  

%%sol   = dde23(    ddefile,    lags,        history,           tspan);
  sol   = dde23(       @dde,    [1],        @history,           tspan);

[t,x_nd] = ode45(@dde_nd,tspan,IC);
  
%% Plots of results  
ax1 = subplot(3,1,1);
        plot(sol.x,sol.y);
    grid on;
    title('sol.x, sol.y')
    xlabel('Time t');    ylabel('x(t)');
 
ax2 = subplot(3,1,2);
        xplot = 1:1:40;
        yplot = deval(sol,xplot,2*n);
    grid on;
    title('deval')
    plot(xplot,yplot);
    xlabel('Time t');    ylabel('x(t)');

ax2 = subplot(3,1,3);
    plot(t,x_nd(:,end)); 
    xlabel('Time t no delay'); ylabel('x(t) no delay');   
    
function xdot=dde(t,x,Z)                % x - state of system (position of each block)
%% INPUT for BIT
z_b=5;                                  %Number of blades
a=0.4445/2;                             %Radius of section#3 in [m]
                                        %Orientation of cutting face =1 at 1st iteration
se=32000;                               %MSE Intrinsic Energy SHOULD BE 32000MPa NEED TO CHECK
e=z_b*a*se;                             %Cutting Energy
%% INPUT for DRILLSTRING
v0=100;                                 %constant speed of drilling                               
n=3;                                    %number of elements for section #3 - length 1350m       
%Mass matrix
m=33299/n;                              %mass of each block
M=m*eye(n);                             %Matrix of drillstring mass
%Stiffness matrix
k=528899/n;                             %spring of each block
K=full(gallery('tridiag',n,-1,2,-1));   %spring
K(end,end)=1;
K=k*K;
K(end,end)=K(end,end)+e;                %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
%Eigenvalue decomposition
[V,D] = eig(K,M);    
%% MVF1                                 
mv=inv(M)*inv(V);
mvzeros=zeros(n,n);
mvzeros(1,1)=mv(1,1);
F1=zeros(n,1); 
F1(1)=v0*t;
MVF1=mvzeros*F1;
%% C damping coefficient
C=534.92;                               %damping
c=C/n;                                  %damping
C=c*eye(n);
%% D
D=D*eye(n);
%% MExlag                               Energy for cutting - Wc only for last element of drillstring (n)
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
        xdot = [x1dot; x2dot]
        
end

function s = history(t)
n=3;
    s=zeros(2*n,1);
    for i=1:2*n
        s(i,1)=1;
    end
end

function xdot_nd=dde_nd(t,x_nd)                % x - state of system (position of each block)
% %% INPUT for BIT
% z_b=5;                                  %Number of blades
% a=0.4445/2;                             %Radius of section#3 in [m]
%                                         %Orientation of cutting face =1 at 1st iteration
% se=32000;                               %MSE Intrinsic Energy SHOULD BE 32000MPa NEED TO CHECK
% e=z_b*a*se;                             %Cutting Energy



% %% INPUT for DRILLSTRING
% v0=100;                                 %constant speed of drilling                               
v0=1;                                   %constant speed of drilling
n=3;                                    %number of elements for section #3 - length 1350m       
% %Mass matrix
% m=33299/n;                              %mass of each block
m=1;                                    %mass of each block
M=m*eye(n);                             %Matrix of drillstring mass
% %Stiffness matrix
% k=528899/n;                             %spring of each block
% K=full(gallery('tridiag',n,-1,2,-1));   %spring
% K(end,end)=1;
% K=k*K;
% K(end,end)=K(end,end)+e;                %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
K=full(gallery('tridiag',n,-1,2,-1));   %spring
K(end,end)=1;                           %put last right bottom value to 1
k=1;                                    %spring of each block
K=k*K;                                  %Stiffness matrix
% %Eigenvalue decomposition
[V,D] = eig(K,M);    
% %% MVF1                                 
f1=v0*t;                                % external force f1=f(t)*M^-1 - 9.81
F1=zeros(n,1);
F1(1)=f1
MVF1=inv(M)*inv(V)*F1
% mv=inv(M)*inv(V);
% mvzeros=zeros(n,n);
% mvzeros(1,1)=mv(1,1);
% F1=zeros(n,1); 
% F1(1)=v0*t;
% MVF1_nd=mvzeros*F1;
% %% C damping coefficient
c=0.2;                                  %damping coefficient
cXdot=zeros(2*n,2*n);
cXdot(1:n, n+1:end)=eye(n);
c=c*eye(n)
% C_nd=534.92;                               %damping
% c_nd=C_nd/n;                                  %damping
% C_nd=c_nd*eye(n);
% %% D
% D_nd=D*eye(n);
DX=zeros(2*n,2*n);
D=D*eye(n)
% %% MExlag                               Energy for cutting - Wc only for last element of drillstring (n)
% E=zeros(n,n);
% E(end,end)=e;                           %last element of vector E is cutting constant e                                       
% ME_nd=E*inv(M);
% %% The differential equations:
%         x1      = zeros(n,1);
%         x2      = zeros(n,1);
%         x1dot   = zeros(n,1);
%         x2dot   = zeros(n,1);
%         xdot    = zeros(2*n,1);
%         x1 = x(1:n);
%         x2 = x(n+1:end);
%         x1dot = x2;
%         x2dot = MVF1-D*x1-c*x2;
%         xdot = [x1dot; x2dot];
    x1_nd      = zeros(n,1);
    x2_nd      = zeros(n,1);
    x1dot_nd   = zeros(n,1);
    x2dot_nd   = zeros(n,1);
    xdot_nd    = zeros(2*n,1);
    x1_nd = x_nd(1:n);
    x2_nd = x_nd(n+1:end);
    x1dot_nd = x2_nd;
    x2dot_nd = MVF1-D*x1_nd-c*x2_nd;
%         x2dot_nd = MVF1_nd - C_nd*x2_nd -D_nd*x1_nd + ME_nd*x1_nd;
         xdot_nd = [x1dot_nd; x2dot_nd];
% 
 end
