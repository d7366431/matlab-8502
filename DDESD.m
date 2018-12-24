clc
clear all

tspan = 1:1:50; %[1 5]
%returns [1,20]
d=[1]; % this is where we specify the vector of taus

% sol   = ddesd(    ddefun,     delays,     history,        tspan);
sol     = dde23(    @ddefun,    d,          @history,       tspan);

y1 = sol.y(1,:); % note the y-solution is a row-wise matrix.
y2 = sol.y(2,:);

% tn  = linspace(1,50,100)
% xn  = deval(sol,tn)

plot(y1,y2)

% sol = 
%      solver: 'dde23'
%     history: @dde_example1/history
%     discont: [0 1 2 3]
%           x: [1x127 double]
%           y: [2x127 double]
%       stats: [1x1 struct]
%          yp: [2x127 double]

% we define the function for the delay. the Y variable is the same as you
% should be used to from an ordinary differential equation. Z is the values
% of all the delayed variables.

    function  dYdt = ddefun(t,Y,Z)
        a =0.25;

        y1 = Y(1);
        y2 = Y(2);

        % Z(:,1) = [y1(t - tau_1); y2(t - tau_2)]
        y1_tau1 = Z(1,1)
        y2_tau1 = Z(2,1)

        dy1dt = y2;
        dy2dt = y2 + y1_tau1;

        dYdt = [dy1dt; dy2dt]
    end

% history function is a column vector giving the value of y1
% and y2 in the past. Here we make these constant values.
% returns [1; 1]
 function y = history(t)
     y = [80; 30];
 end