%% ODE to solve y'(t) = y(t)' - 2y(t - 1)
clc
clear all

lags = 1; 
tspan = 1:1:100

%%sol   = dde23(    ddefile,    lags,       history,            tspan);
  sol   = dde23(    @dde,       lags,       @history,           tspan)

xais    = linspace(1,100,1000)
yaxis   = deval(sol,xais)

plot(xais,yaxis)


function dydt = dde(t,y,s)
dydt = y - 2*s;
end

function s = history(t)
s = t - 1;
end

% Syntax sxint = deval(sol,xint)
% sxint = deval(sol,xint) evaluates the solution of a differential equation problem at each element of the
% vector xint. For each i, sxint(:,i) is the solution corresponding to xint(i).
% The input argument sol is a structure returned by an initial value problem solver (ode45, ode23, ode113,
% ode15s, ode23s, ode23t, ode23tb) or the boundary value problem solver (bvp4c). The ordered row vector
% sol.x contains the independent variable. For each i, the column sol.y(:,i) contains the solution at sol.x(i).
% The structure sol also contains data needed for interpolating the solution in (x(i),x(i+1)). The form of this
% data depends on the solver that creates sol. The field sol.solver contains the name of that solver.
% Elements of xint must be in the interval [sol.x(1),sol.x(end)].