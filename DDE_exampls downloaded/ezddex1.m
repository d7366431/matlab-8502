function ezddex1
%EZDDEX1  Example 1 for DDE23 solved with different syntax using EZDDE23.
%   This is a simple example of Wille' and Baker that illustrates the
%   straightforward formulation, computation, and plotting of the solution 
%   of a system of delay differential equations (DDEs). 
%
%   The differential equations
%
%        y'_1(t) = y_1(t-1)  
%        y'_2(t) = y_1(t-1)+y_2(t-0.2)
%        y'_3(t) = y_2(t)
%
%   are solved on [0, 5] with history y_1(t) = 1, y_2(t) = 1, y_3(t) = 1 
%   for t <= 0. 
%
%   The lags are specified as a vector [1, 0.2]. Generally the history is
%   evaluated by a function, but when the history is constant, as it is 
%   here, it can be supplied as a vector.  EZDDE23 differs from DDE23 only
%   in the way that the delay differential equations are coded.  As seen in
%   the subfunction DDEs, terms with delays are represented by arguments
%   YLAGJ which correspond to the column vectors y(T - LAGS(J)).  
%
%   See also DDE23, FUNCTION_HANDLE.

history = @(t) ones(3,1);  % Or: history = [1;1;1];

sol = ezdde23(@DDEs,[1, 0.2],history,[0, 5]);
figure;
plot(sol.x,sol.y)
title('An example of Wille'' and Baker.'); 
xlabel('time t');
ylabel('solution y');

% -------------------------------------------------------------------------

function dydt = DDEs(t,y,ylag1,ylag2)
% Differential equations function for EZDDEX2.
dydt = [ ylag1(1)
         ylag1(1) + ylag2(2)
         y(2)               ];
