n=2;                                    %number of elements for section #3 - length 1350m     
IC = zeros(2*n,1)
tspan = 0:0.001:1;
[t,x] = ode45(@fr,tspan,IC);

plot(t,x(:,end)); 
xlabel('t'); ylabel('x'); 

%% Create function that returns sdot input is time (t) and state state (x)
function xdot=fr(t,x)                   % x - state of system (position of each block)
%% INPUT for BIT
z=5;                                    %number of blades
a=0.4445/2;                             %radius of section#3 in [m]
                                        %Orientation of cutting face =1 at 1st iteration
se=1;                                   %MSE Intrinsic Energy SHOULD BE 20 NEED TO CHECK???????????
e=z*a*se;                               %Cutting Energy
%% INPUT for DRILLSTRING
v0=1;                                   %constant speed of drilling
c=0.2;                                  %damping coefficient
n=2;                                    %number of elements for section #3 - length 1350m       
m=1;                                    %mass of each block
M=m*eye(n)                              %Matrix of drillstring mass
k=1;
K=full(gallery('tridiag',n,-1,2,-1));   %spring
K=k*K;
K(end,end)=k+e;                          %put last right bottom value as sum of spring (K) and "Cutting constant" (e)
                                        %spring of each block
                                        %Stiffness matrix
[V,D] = eig(K,M);                        %Eigenvalue decomposition
%% MVF1
mv=inv(M)*inv(V);
mvzeros=zeros(n,n);
mvzeros(1,1)=mv(1,1);
F1=zeros(n,1); 
F1(1)=v0*t
MVF1=mvzeros*F1
% f1=v0*t;                                % external force f1=f(t)*M^-1 - 9.81
% F1=zeros(n,1);
% F1(1)=f1
% MVF1=inv(M)*inv(V)*F1
%% Dx
DX=zeros(2*n,2*n);
D=D*eye(n)
%% cx
cXdot=zeros(2*n,2*n);
cXdot(1:n, n+1:end)=eye(n);
c=c*eye(n)
%% MEx(t-tau)                           Energy for cutting - Wc only for last element of drillstring (n)
E=zeros(n,n);
E(end,end)=e;                            %last element of vector E is cutting constant e                                       
MEsm=inv(M)*E
%ME=zeros(n,n)
%ME(end, end)=MEsm

%% The differential equations:
%[x]=[x1;x2]=[x;xdot]
%[xdot]=[x2;big equation]
        x1      = zeros(n,1);
        x2      = zeros(n,1);
        x1dot   = zeros(n,1);
        x2dot   = zeros(n,1);
        xdot    = zeros(2*n,1);
        x1 = x(1:n)
        x2 = x(n+1:end,1)
        x1dot = x2;
        x2dot = MVF1-D*x1-c*x2-MEsm*x1;
        xdot = [x1dot; x2dot];
end