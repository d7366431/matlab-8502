function ezddex3
%EZDDEX3  Example 3 for DDESD solved with different syntax using EZDDESD.
%   This example solves a system of two delay differential equations with
%   a state-dependent delay
%       y1'(t) =  y2(t)  
%       y2'(t) = -y2(exp(1-y2(t))) * y2(t)^2 * exp(1-y2(t))  
%   The system has the analytical solution
%       y1(t) = log(t)
%       y2(t) = 1/t
%   which serves as the history for t < 0.1.
%   
%   The Problem comes from W.H. Enright and H. Hayashi, The Evaluation 
%   of Numerical Software for Delay Differential Equations, pp. 179-192 
%   in R. Boisvert (Ed.), The Quality of Numerical Software: Assessment 
%   and Enhancement, Chapman & Hall, London, 1997.
%  
%   See also DDESD, FUNCTION_HANDLE.

t0 = 0.1;
tfinal = 5;
tspan = [t0, tfinal];

sol = ezddesd(@DDEs,@ddex3delay,@ddex3hist,tspan);

% Exact solution
texact = linspace(t0,tfinal);
yexact = ddex3hist(texact);

figure
plot(texact,yexact,sol.x,sol.y,'o')
legend('y_1, exact','y_2, exact','y_1, ezddesd','y_2, ezddesd')
xlabel('time t')
ylabel('solution y')
title('D1 problem of Enright and Hayashi')

% -----------------------------------------------

function v = ddex3hist(t)
% History function for EZDDEX3.  
  v = [ log(t); 1./t];

% -----------------------------------------------

function d = ddex3delay(t,y)
% State dependent delay function for EZDDEX3.  
  d = exp(1 - y(2));

% -----------------------------------------------

function dydt = DDEs(t,y,ylag)
% Differential equations function for EZDDEX3.    
  dydt = [ y(2); -ylag(2)*y(2)^2*exp(1 - y(2))];
