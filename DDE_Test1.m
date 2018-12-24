%% ODE to solve y?(t) = y(t) ? 2y(t ? 1)
clc
clear all
tspan = 1:1:5;
%%sol   = dde23(    ddefile,    lags,       history,            tspan);
  sol   = dde23(       @dde,    [1],        @history,           tspan);
plot(sol.x,sol.y);
xlabel('time t');
ylabel('y(t)');

function xdot = dde(t,x,Z) 
    xlag=Z(:,1)
    n=3;
    x1 = x(1:n);
    x2 = x(n+1:end);
    x1dot = 2*x1 + xlag;
    x2dot = x2;
    xdot = [x1dot; x2dot];
end

function s = history(t)
n=3
    s=zeros(n,1)
    for i=1:n
        s(i,1)=1
    end
    %s=s(i)
end