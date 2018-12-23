clc
clear all

%%  sol = dde23(    ddefile,    lags,       history,    tspan);
    sol = dde23(    @exam1f,    [1, 0.2],   ones(3,1),  [0, 5]);
    
plot(sol.x,sol.y);
title('Figure 1. Example 3 of Wille'' and Baker.')
xlabel('time t');
ylabel('y(t)');

function v = exam1f(t,y,Z)
    ylag1 = Z(:,1)
    ylag2 = Z(:,2)
    v = zeros(3,1)
    v(1) = ylag1(1)
    v(2) = ylag1(1) + ylag2(2)
    v(3) = y(2)
end