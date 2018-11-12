function cyi = f16_cy(beta_deg, ail_deg, rud_deg)
% Data and interpolation for F-16 z-Force Coefficient.
% 
% Inputs:
%   beta_deg    sideslip angle, in degrees
%   ail_deg     aileron deflection, in degrees
%   rud_deg     rudder deflection, in degrees
%
% Outputs:
%   cyi     lateral force coefficient
%
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Limit the aileron deflection to +/- 20 degrees.
ail_ltd_deg = min(20,max(-20,ail_deg));

% Limit the rudder deflection to +/- 30 degrees.
rud_ltd_deg = min(30,max(-30,rud_deg));

% Compute the lateral force coefficient
cyi = -0.02*beta_deg + 0.021*(ail_ltd_deg/20) + 0.086*(rud_ltd_deg/30);