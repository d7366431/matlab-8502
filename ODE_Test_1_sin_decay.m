clear all
clc 

n=3;                                    %number of elements for section #3 - length 1350m     
IC = zeros(2*n,1)
tspan = 0:0.001:50;
[t,x] = ode45(@fr,tspan,IC);

plot(t,x(:,end)); 
xlabel('t'); ylabel('x'); 

%% Create function that returns sdot input is time (t) and state state (x)
function xdot=fr(t,x)                   % x - state of system (position of each block)
% INPUT for DRILLSTRING
v0=10;                                   %constant speed of drilling
c=0.2;                                  %damping coefficient
n=3;                                    %number of elements for section #3 - length 1350m       
m=1;                                    %mass of each block
M=m*eye(n)                              %Matrix of drillstring mass
K=full(gallery('tridiag',n,-1,2,-1));   %spring
K(end,end)=1;                           %put last right bottom value to 1
k=1;                                    %spring of each block
K=k*K;                                  %Stiffness matrix
[V,D] = eig(K,M)                        %Eigenvalue decomposition
%% MVF1
f1=v0*t;                                % external force f1=f(t)*M^-1 - 9.81
F1=zeros(n,1);
F1(1)=f1
MVF1=inv(M)*inv(V)*F1
%% Dx
DX=zeros(2*n,2*n);
D=D*eye(n)
%% cx
cXdot=zeros(2*n,2*n);
cXdot(1:n, n+1:end)=eye(n);
c=c*eye(n)
%% The differential equations:
%[x]=[x1;x2]=[x;xdot]
%[xdot]=[x2;big equation]
        x1      = zeros(n,1);
        x2      = zeros(n,1);
        x1dot   = zeros(n,1);
        x2dot   = zeros(n,1);
        xdot    = zeros(2*n,1);
        x1 = x(1:n);
        x2 = x(n+1:end);
        x1dot = x2;
        x2dot = MVF1-D*x1-c*x2;
        xdot = [x1dot; x2dot];
end