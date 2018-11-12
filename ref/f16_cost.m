function [cost, xd, x, u] = f16_cost(ux0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the cost used in assessing a vehicle state as
% being in trim.
%
%
% Inputs:
%   ux0     [5x1] vector, initial conditions on the inputs and states
%       ux0(1): Beta (rad)
%       ux0(2): Elevator deflection (deg)
%       ux0(3): Aileron deflection (deg)
%       ux0(4): Rudder deflection (deg)
%       ux0(5): Throttle (0-1)
%
%
% Outputs:
%   cost    scalar, converged cost of the trimmed system
%   xd      n-element vector, derivative of state vector
%   x       n-element vector, trimmed state vector
%   u       m-element vector, trimmed control vector
%
%
% Version history:
%   2011.09.08  KRB  Initial release
%   2016.10.05  KRB  Use alpha as a fixed flight condition parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Vt_fps alt_ft alpha_rad

% Compute the power.
pow = f16_tgear(ux0(5));

% Build the state : Vt, AOA, AOS, Phi, Theta, Psi, P, Q, R, NED, pow.
x = [Vt_fps;
     alpha_rad;
     ux0(1);
     zeros(8,1);
     alt_ft;
     pow];


% Apply any designated constraints.
%x = f16_constr(x,trim);

% Build the control vector : thtl, de, da, dr.
u = [ux0(5);ux0(2);ux0(4);ux0(3)];


% Compute the state derivatives.
xd = f16_nonlinear_model(0, x, u);


% Compute the cost of this vector.  Desire this cost to be small.  The cost
% is a weighted sum of the pertinent states squared.
weight = [1;100;100;10;10;10];

cost = weight'*[xd(1);xd(2);xd(3);xd(7);xd(8);xd(9)].^2;
