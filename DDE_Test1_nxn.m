clc
clear all

tspan = 1:1:5;
%%sol   = dde23(    ddefile,    lags,        history,           tspan);
  sol   = dde23(       @dde,    [1],        @history,           tspan);

plot(sol.x,sol.y);
xlabel('time t');
ylabel('x(t)');

function qdot=dde(t,x,Z)                % x - state of system (position of each block)
%% INPUT for BIT
z=5;                                    %number of blades
a=0.4445/2;                             %radius of section#3 in [m]
                                        %Orientation of cutting face =1 at 1st iteration
se=1;                                   %MSE Intrinsic Energy SHOULD BE 20 NEED TO CHECK???????????
e=z*a*se;                               %Cutting Energy

%% INPUT for DRILLSTRING
v0=10;                                   %constant speed of drilling
c=0.2;                                  %damping coefficient
n=3;                                    %number of elements for section #3 - length 1350m       
m=1;                                    %mass of each block
M=m*eye(n);                             %Matrix of drillstring mass
k=1;
K=full(gallery('tridiag',n,-1,2,-1));   %spring
K=k*K;
K(end,end)=k+e;                         %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
                                        %spring of each block
                                        %Stiffness matrix
[V,D] = eig(K,M);                       %Eigenvalue decomposition
%% MVF1
f1=v0*t; 
F1=zeros(n,1); 
F1(1)=f1;
mv=inv(M)*inv(V);
mvF1=mv*F1;
MVF1=[zeros(n,1);mvF1]
%% Dx
Dn=zeros(2*n,2*n);
Dn(n+1:end, 1:n)=D*eye(n)
%% cxdot
C=zeros(2*n,2*n);
C(1:n, n+1:end)=eye(n);
C(n+1:end, n+1:end)=c*eye(n)
%% MEx(t-tau)                           Energy for cutting - Wc only for last element of drillstring (n)
E=zeros(n,n);
E(end,end)=e;                           %last element of vector E is cutting constant e                                       
MEsm=inv(M)*E;
ME=zeros(2*n,2*n);
ME(n+1:end, 1:n)=MEsm
%% The differential equations:
%[x]    =   [x1;x2]=[x;xdot]
%[xdot] =   [x2;big equation]

xdot = MVF1 - Dn*x + ME*x - 1%C*xlag;
q = zeros(2*n,2*n);
q(n+1:end,n+1:end)=inv(V);
qdot=q*xdot;

end

function s = history(t)
n=3
    s=zeros(n,1)
    for i=1:n
        s(i,1)=1
    end
    %s=s(i)
end