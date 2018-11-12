%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AENG - 555 : Guidance and Control of Aerospace Vehicles
% Saint Louis University
% Ken Buckholtz, Ph.D.
%
% Version history:
%   2011.09.08  KRB  Initial Release
%   2016.10.05  KRB  Use alpha_rad as a fijxed flight condition for trim
%
% This script trims the nonlinear F-16 model at the prescribed flight
% condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
format short g
global Vt_fps alt_ft xcg alpha_rad

% Flight condition.
Vt_fps = 500; 
alt_ft = 1;
alpha_rad = 2*(pi/180);
    
% Define xcg location.
xcg = 0.30;

% Initial "guess" for trim routine.
beta_rad  = 0;
de_deg    = -1.931;
dr_deg    = 0;
da_deg    = 0;
thtl      = 0.1485;

% Initial guess for free parameters.
ux0 = [beta_rad; de_deg; da_deg; dr_deg; thtl];

% Initializing optimization options and running optimization:
OPTIONS = optimset('TolFun',1e-10, ...
                   'TolX',1e-10,...
                   'MaxFunEvals',5e+04,...
                   'MaxIter',1e+04);

% Perform the optimization to find the vector UX that minimizes the 
% cost.
[UX,FVAL,EXITFLAG,OUTPUT] = fminsearch('f16_cost',ux0,OPTIONS);

if EXITFLAG
    disp('Trim algorithm converged.');
else
    disp('Trim algorithm DID NOT converge.');
end

% Evaluate the const using the optimized vector UX to get the cost and
% associated state vector, x, and control vector, u, that trims the 
% model.
[cost, xd, x, u] = f16_cost(UX);
    
