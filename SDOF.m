
%% INPUT bottom of section 3 - 444.5mm
m=40363         %mass of drillstring
t=0.001         %time step
T=20           %end modeling time
c=2715          %damping
k=528889        %spring


%% SYSTEM PARAMETERS
w0=sqrt(k/m);                            %natural frquency of system
theta=c/(2*sqrt(k*m));                   %decay coefficient
wd=sqrt(pow2(w0)-pow2(theta*w0));        %wd of system

%% TIME STEPS
T_vector = (-T:t:T);                    %from -T untill T with time step t
Q_function = zeros(1,length(T_vector));
S1 = zeros(1,length(T_vector));

%% LOOP
S1_position = zeros(1,length(T_vector));
S1_velocity = zeros(1,length(T_vector));
S1_acceleration = zeros(1,length(T_vector));

for i=1:length(T_vector)
    if  T_vector(i)<0;
        S1_position(i) =0;
        S1_velocity(i) =0;
        S1_acceleration(i)=0;
        Q_function(i) = 0;
    else       
        S1_position(i) =(exp(-theta*w0*T_vector(i)))*(sin(wd*T_vector(i)))*(w0^2/wd);
        S1_vel_term1(i) = wd*(exp(-theta*w0*T_vector(i)))*cos(wd*T_vector(i));
        S1_vel_term2(i) = theta*w0*(exp(-theta*w0*T_vector(i)))*sin(wd*T_vector(i));
        S1_velocity(i) = (pow2(w0)/wd)*(S1_vel_term1(i)-S1_vel_term2(i));
        S1_accel_term1(i)=theta*w0*exp(-theta*w0*T_vector(i));
        S1_accel_term2(i)=(theta*w0*sin(wd*T_vector(i)))-wd*cos(wd*T_vector(i));
        S1_accel_term3(i)=S1_accel_term1(i)*S1_accel_term2(i);
        S1_accel_term4(i)= -wd*exp(-theta*w0*T_vector(i));
        S1_accel_term5(i)=wd*sin(wd*T_vector(i))+theta*cos(wd*T_vector(i));
        S1_accel_term6(i)=S1_accel_term4(i)*S1_accel_term5(i);
        S1_acceleration(i)=(S1_accel_term3(i)+S1_accel_term6(i))*(pow2(w0)/wd);
        Q_function(i) = 3*sin((2*3.14/10)*T_vector(i)); %sin excitation force
        %Q_function(i)=9.81*m*(pow2(T_vector(i)))/2;
    end
end

P=conv(S1_position,Q_function, 'same')*t;
V=conv(S1_velocity,Q_function, 'same')*t;
A=conv(S1_acceleration,Q_function, 'same')*t;

%% PLOTS
ax1 = subplot(5,1,1);
plot(T_vector,Q_function) %,T_vector,P
grid on;
title('Q-function, Position')
xlim([0,inf])

ax2 = subplot(5,1,2);
plot(T_vector,S1_velocity)
grid on;
title('S1-velocity')
xlim([0,inf])

ax3 = subplot(5,1,3);
plot(T_vector,P,T_vector,Q_function)
grid on;
title('Position')
xlim([0,inf])
legend('P','Q function');

ax2 = subplot(5,1,4);
plot(T_vector,V,T_vector,Q_function);
grid on;
title('Velocity, Acceleration');
xlim([0,inf])
legend('V','Q function');

ax5 = subplot(5,1,5);
plot(T_vector,A,T_vector,V);
grid on;
title('Acceleration');
xlim([0,inf])
legend('A','V');
%linkaxes([ax1,ax2,ax3,ax4,ax5],'xy')