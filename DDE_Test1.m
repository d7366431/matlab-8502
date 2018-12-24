%% ODE to solve y?(t) = y(t) ? 2y(t ? 1)
clc
clear all

tspan = 1:1:10;

%%sol   = dde23(    ddefile,    lags,       history,            tspan);
  sol   = dde23(       @dde,    [1],       [1; 1; 1],           [1 5]);
%@history is a function that returns a solution for the system at t<=t0, specified in one of these ways:
            % -> A function of t such that y = history(t) returns the solution y(t) for t ? t0 as a column vector
            % -> A constant column vector, if y(t) is constant
            % -> The solution sol from a previous integration, if this call continues that integration <tspan> 
            %    specifies t0 and tend for your solution.

% xais    = linspace(1,5,100);
% yaxis   = deval(sol,xais);
% 
% plot(xais,yaxis)

plot(sol.x,sol.y);
xlabel('time t');
ylabel('y(t)');
%%************************************************************************

%@ddefile is ODE function handle <lags> is an array of constants specifying the delay for each variable in your function
function xdot = dde(t,x,Z) 
xlag=Z(:,1)
n=3;
x1 = x(1:n);
x2 = x(n+1:end);
x1dot = 2*x1 + xlag
x2dot = x2
xdot = [x1dot; x2dot];
end

function s = history(t)
s = [1;1];
end