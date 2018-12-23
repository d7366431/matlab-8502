clear all
clc

%% INPUT
n=1000;                                 %number of elements for section #3 - length 1350m

M=33299;                                %mass of drillstring
m=M/n;                                  %mass of each block

C=534.92;                               %damping
c=C/n;                                  %damping

K=528899;                               %spring
k= n*K;                                 %Spring coefficient

%% TIME
t=0.0001;                               %time step
T=20;                                   %end modeling time
T_vector = -T:t:T

%% SYSTEM PARAMETERS
M=m*eye(n)                              %Matrix of drillstring mass

K=full(gallery('tridiag',n,-1,2,-1));   %spring
K(end,end)=1  ;                         %put last right bottom value to 1
K=k*K  ;                                %Stiffness matrix

                                        %Eigenvalue decomposition                                        
[V,D] = eig(K,M);                       %produces a diag matrix D of generalized eigenvalues 
                                        %and a full matrix V whose columns are the corresponding 
                                        %eigenvectors so that A*V = B*V*D. 
                                        %The eigenvectors are scaled so that the norm of each is 1.0.  
                                        
w0=D.^0.5;                              
w0=diag(w0) ;                           %Natural frequencies (sqrtDii)

theta=c./(2*m*w0);                      %damping coefficient

wd=(sqrt(1-theta.^2)).*w0;              %wd of system

%% DIRECT SOLUTION FOR STEP FUNCTION
%sum_d = 0;
T_vector_2=0:0.0001:T;
x_d = zeros(length(T_vector_2),1);
n2=500;
j=n2;
Qs=1;
A1=1;
%Bottom of the string
value=n-1;

for z=1:length(T_vector_2)
    sum_d=0;
    for i=1:n2
        sum_d=sum_d+((V(n2,i)*V(1,i))/((w0(i)^2)*sqrt(1-theta(i)^2)))*sin(wd(i)*T_vector_2(z)+acos(theta(i)));
    end
    x_d(z,1) =Qs-Qs*0.9056*(n2^2)*exp((-c/2)*T_vector_2(z))*sum_d;
end

%% SOLUTION FOR CONVOLUTION
sum = 0;
X = zeros(length(T_vector),1);
Q = zeros(length(T_vector),1);
j=n;
Qs=1;

for u = 1:length(T_vector)
    if  T_vector(u) < 0
        X(u,1) = zeros(1,1);
        Q(u) = 0;
    else
        for i = 1:n
            sum = sum + (k*(V(j,i)*V(1,i)) /  ((w0(i))*sqrt(1-theta(i)^2))) * sin(wd(i)*T_vector(u));
        end
        s = exp( -(c/(2*m))*T_vector(u))*sum;
        sum = 0;
        X(u,1) = s;
        Q_sinwave(u)=3*sin((2*3.14/10)*T_vector(u));
        %Q_sinwave(u) = 3*sin(0.628*T_vector(u));
        Q_gravity(u)=(9.81*(m)*(T_vector(u)^2))/2;
    end
end

%%
Y_sinwave = conv(Q_sinwave,X,'same')*t;
Y_gravity = conv(Q_gravity,X,'same')*t;

%% PLOTS
ax0 = subplot(3,1,1);
plot(T_vector_2,x_d)
grid on;
title('Direct Solution for Step Function')
xlim([0,inf])
ylabel('Y');
xlabel('time [s]');

ax1 = subplot(3,1,2);
plot(T_vector,Q_sinwave,T_vector,Y_sinwave)
grid on;
title('Q sinwave,Y sinwave')
xlim([0,inf])
ylabel('Y');
xlabel('time [s]');

ax1 = subplot(3,1,3);
plot(T_vector,Q_gravity,T_vector,Y_gravity)
grid on;
title('Q gravity, Y gravity')
xlim([0,inf])
ylabel('Y');
xlabel('time [s]');