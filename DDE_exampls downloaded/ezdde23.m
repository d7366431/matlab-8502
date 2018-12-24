function sol = ezdde23(ddefun,lags,history,tspan,options,varargin) 
%EZDDE23  Solve delay differential equations (DDEs) with constant delays.
%   SOL = EZDDE23(DDEFUN,LAGS,HISTORY,TSPAN) integrates a system of DDEs 
%   y'(t) = f(t,y(t),y(t - tau_1),...,y(t - tau_k)). The constant, positive 
%   delays tau_1,...,tau_k are input as the vector LAGS. DDEFUN is a function 
%   handle. DDEFUN(T,Y,YLAG1,YLAG2,...,YLAGK) must return a column vector 
%   corresponding to f(t,y(t),y(t - tau_1),...,y(t - tau_k)). In the call to 
%   DDEFUN, a scalar T is the current t, a column vector Y approximates y(t), 
%   and a column vector YLAGJ approximates y(t - tau_j) for delay 
%   tau_j = LAGS(J). The DDEs are integrated from T0 to TF where T0 < TF and 
%   TSPAN = [T0 TF]. The solution at t <= T0 is specified by HISTORY in one 
%   of three ways: HISTORY can be a function handle, where for a scalar T, 
%   HISTORY(T) returns a column vector y(t). If y(t) is constant, HISTORY 
%   can be this column vector. If this call to EZDDE23 continues a previous 
%   integration to T0, HISTORY can be the solution SOL from that call.
%
%   EZDDE23 produces a solution that is continuous on [T0,TF]. The solution is
%   evaluated at points TINT using the output SOL of EZDDE23 and the function
%   DEVAL: YINT = DEVAL(SOL,TINT). The output SOL is a structure with 
%       SOL.x  -- mesh selected by EZDDE23
%       SOL.y  -- approximation to y(t) at the mesh points of SOL.x
%       SOL.yp -- approximation to y'(t) at the mesh points of SOL.x
%       SOL.solver -- 'dde23'  (EZDDE23 is a driver for DDE23.)
%
%   SOL = EZDDE23(DDEFUN,LAGS,HISTORY,TSPAN,OPTIONS) solves as above with default
%   parameters replaced by values in OPTIONS, a structure created with the
%   DDESET function. See DDESET for details. Commonly used options are
%   scalar relative error tolerance 'RelTol' (1e-3 by default) and vector of
%   absolute error tolerances 'AbsTol' (all components 1e-6 by default).
%
%   EZDDE23 can solve problems with discontinuities in the solution prior to T0
%   (the history) or discontinuities in coefficients of the equations at known
%   values of t after T0 if the locations of these discontinuities are
%   provided in a vector as the value of the 'Jumps' option.
%
%   By default the initial value of the solution is the value returned by
%   HISTORY at T0. A different initial value can be supplied as the value of
%   the 'InitialY' property. 
%
%   With the 'Events' property in OPTIONS set to a function handle EVENTS, 
%   EZDDE23 solves as above while also finding where event functions 
%   g(t,y(t),y(t - tau_1),...,y(t - tau_k)) are zero. For each function 
%   you specify whether the integration is to terminate at a zero and whether 
%   the direction of the zero crossing matters. These are the three column 
%   vectors returned by EVENTS: 
%     [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y,YLAG1,...,YLAGK). 
%   For the I-th event function: VALUE(I) is the value of the function, 
%   ISTERMINAL(I) = 1 if the integration is to terminate at a zero of this 
%   event function and 0 otherwise. DIRECTION(I) = 0 if all zeros are to
%   be computed (the default), +1 if only zeros where the event function is
%   increasing, and -1 if only zeros where the event function is decreasing. 
%   The field SOL.xe is a row vector of times at which events occur. Columns
%   of SOL.ye are the corresponding solutions, and indices in vector SOL.ie
%   specify which event occurred.   
%   
%   Examples
%     EZDDE23 is used exactly like DDE23 except for the way the DDEs are
%     coded.  EZDDEX1 and EZDDEX2 are the DDEX1 and DDEX2 examples of 
%     DDE23 coded with the alternative syntax.
%
%   Class support for inputs TSPAN, LAGS, HISTORY, and the result of 
%     DDEFUN(T,Y,YLAG1,...,YLAGK): 
%     float: double, single
%
%   See also DDESET, DDEGET, DEVAL.

%   EZDDE23 tracks discontinuities and integrates with the explicit Runge-Kutta
%   (2,3) pair and interpolant of ODE23. It uses iteration to take steps
%   longer than the lags.

%   Details are to be found in Solving DDEs in MATLAB, L.F. Shampine and
%   S. Thompson, Applied Numerical Mathematics, 37 (2001). 

% Get number of equations and number of lags.
if isnumeric(history)
  y0 = history;
elseif isstruct(history)
  y0 = history.y(:,end);
else
  y0 = feval(history,tspan(1),varargin{:});
end 
neq = length(y0);
nlags = length(lags);

% Initialize.
if nlags ~= 0
    onevec = ones(1,nlags);
else
    ZC = cell(1);
end
if nargin < 5, options = []; end

userEventFcn = ddeget(options,'Events',[]);
if ~isempty(userEventFcn)
    options = ddeset(options,'Events',@wrappedEventFcn);
end

sol = dde23(@ezdde,lags,history,tspan,options,varargin{:});


%================================================

function dydt = ezdde(t,y,Z,varargin)
    if nlags ~= 0
        ZC = mat2cell(Z,neq,onevec);
    end
    dydt = feval(ddefun,t,y,ZC{:},varargin{:});
end % ezdde

%------------------------------------------------
function [value,isterminal,direction] = wrappedEventFcn(t,y,Z,varargin)
    if nlags ~= 0
        ZC = mat2cell(Z,neq,onevec);
    end
    [value,isterminal,direction] = feval(userEventFcn,t,y,ZC{:},varargin{:});
end % wrappedEventFcn

%================================================

end % ezdde23