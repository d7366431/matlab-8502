function sol = ezddesd(ddefun,delays,history,tspan,options) 
%EZDDESD  Solve delay differential equations (DDEs) with general delays.
%   SOL = EZDDESD(DDEFUN,DELAYS,HISTORY,TSPAN) integrates a system of DDEs 
%   y'(t) = f(t,y(t),y(d(1)),...,y(d(k))). The delays d(j) can depend on 
%   both t and y(t).  DDEFUN and DELAYS are function handles. DELAYS(T,Y)  
%   must return a column vector of delays d(j). EZDDESD imposes the requirement 
%   that d(j) <= t by using min(d(j),t).  DDEFUN(T,Y,YLAG1,YLAG2,...,YLAGK)
%   must return a column vector corresponding to f(t,y(t),y(d(1)),...,y(d(k))).
%   In the call to DDEFUN, a scalar T is the current t, a column vector Y 
%   approximates y(t), and a column vector YLAGJ approximates y(d(j)) for 
%   delay d(j). The DDEs are integrated from T0 to TF where T0 < TF and 
%   TSPAN = [T0 TF]. The solution at t <= T0 is specified by HISTORY in one 
%   of three ways: HISTORY can be a function handle, where for a scalar T, 
%   HISTORY(T) returns the column vector y(t). If y(t) is constant, HISTORY 
%   can be this column vector. If this call to EZDDESD continues a previous 
%   integration to T0, HISTORY can be the solution SOL from that call.
%
%   EZDDESD produces a solution that is continuous on [T0,TF]. The solution is
%   evaluated at points TINT using the output SOL of EZDDESD and the function
%   DEVAL: YINT = DEVAL(SOL,TINT). The output SOL is a structure with 
%       SOL.x  -- mesh selected by EZDDESD
%       SOL.y  -- approximation to y(t) at the mesh points of SOL.x
%       SOL.yp -- approximation to y'(t) at the mesh points of SOL.x
%       SOL.solver -- 'ddesd' (EZDDESD is a driver for DDESD.)
%
%   SOL = EZDDESD(DDEFUN,DELAYS,HISTORY,TSPAN,OPTIONS) solves as above with 
%   default parameters replaced by values in OPTIONS, a structure created 
%   with the DDESET function. See DDESET for details. Commonly used options 
%   are scalar relative error tolerance 'RelTol' (1e-3 by default) and vector
%   of absolute error tolerances 'AbsTol' (all components 1e-6 by default).
%
%   By default the initial value of the solution is the value returned by
%   HISTORY at T0. A different initial value can be supplied as the value of
%   the 'InitialY' property. 
%
%   With the 'Events' property in OPTIONS set to a function handle EVENTS, 
%   EZDDESD solves as above while also finding where event functions 
%   g(t,y(t),y(d(1)),...,y(d(k))) are zero. For each function you specify 
%   whether the integration is to terminate at a zero and whether the 
%   direction of the zero crossing matters. These are the three vectors 
%   returned by EVENTS: 
%     [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y,YLAG1,YLAG2,...,YLAGK). 
%   For the I-th event function: VALUE(I) is the value of the function, 
%   ISTERMINAL(I) = 1 if the integration is to terminate at a zero of this 
%   event function and 0 otherwise. DIRECTION(I) = 0 if all zeros are to
%   be computed (the default), +1 if only zeros where the event function is
%   increasing, and -1 if only zeros where the event function is decreasing. 
%   The field SOL.xe is a row vector of times at which events occur. Columns
%   of SOL.ye are the corresponding solutions, and indices in vector SOL.ie
%   specify which event occurred.   
%
%   If all the delay functions have the form d(j) = t - tau_j, you can set 
%   the argument DELAYS to a constant vector DELAYS(j) = tau_j. With delay 
%   functions of this form, EZDDESD is used exactly like EZDDE23.
%   
%   Example
%     EZDDESD is used exactly like DDESD except for the way that the DDEs 
%     are coded.  EZDDEX3 is the DDEX3 example of DDESD coded with the 
%     alternative syntax.
%
%   Class support for inputs TSPAN, HISTORY, and the results of DELAYS(T,Y) 
%   and DDEFUN(T,Y,YLAG1,YLAG2,...,YLAGK):
%     float: double, single
%
%   See also EZDDE23, DDESET, DDEGET, DEVAL.

% Get number of equations and number of delays.
if isnumeric(history)
  y0 = history;
elseif isstruct(history)
  y0 = history.y(:,end);
else
  y0 = feval(history,tspan(1));
end 
neq = length(y0);
if isnumeric(delays)
    ndelays = length(delays);
else
    ndelays = length(feval(delays,tspan(1),y0));
end

% Initialize.
if ndelays ~= 0
    onevec = ones(1,ndelays);
else
    ZC = cell(1);
end
if nargin < 5, options = []; end
userEventFcn = ddeget(options,'Events',[]);
if ~isempty(userEventFcn)
    options = ddeset(options,'Events',@wrappedEventFcn);
end

sol = ddesd(@ezdde,delays,history,tspan,options);

%================================================
function dydt = ezdde(t,y,Z)
    if ndelays ~= 0
        ZC = mat2cell(Z,neq,onevec);
    end
    dydt = feval(ddefun,t,y,ZC{:});
end % ezdde

%------------------------------------------------
function [value,isterminal,direction] = wrappedEventFcn(t,y,Z,varargin)
    if ndelays ~= 0
        ZC = mat2cell(Z,neq,onevec);
    end
    [value,isterminal,direction] = feval(userEventFcn,t,y,ZC{:},varargin{:});
end % wrappedEventFcn

%================================================

end % ezddesd
