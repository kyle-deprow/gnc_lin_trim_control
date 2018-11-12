function cli = f16_cl(alpha_deg, beta_deg)
% Data and interpolation for F-16 Rolling Moment Coefficient.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip, in degrees
%
% Outputs:
%   cli     rolling moment coefficient
% 
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Independent variables' breakpoints.
alphtab_deg = -10:5:45;
betatab_deg = 0:5:30;

% Bound aoa and aos to the table ranges.
alpha_ltd_deg = max(-10, min(45, alpha_deg));
beta_ltd_deg  = max(-30, min(30, beta_deg));

% Cl table
cltabl = [ zeros(1,12); ...
          -0.001, -0.004, -0.008, -0.012, -0.016, -0.019, -0.020, -0.020, ...
          -0.015, -0.008, -0.013, -0.015;
          -0.003, -0.009, -0.017, -0.024, -0.030, -0.034, -0.040, -0.037, ...
          -0.016, -0.002, -0.010, -0.019;
          -0.001, -0.010, -0.020, -0.030, -0.039, -0.044, -0.050, -0.049, ...
          -0.023, -0.006, -0.014, -0.027;
           0.000, -0.010, -0.022, -0.034, -0.047, -0.046, -0.059, -0.061, ...
          -0.033, -0.036, -0.035, -0.035;
           0.007, -0.010, -0.023, -0.034, -0.049, -0.046, -0.068, -0.071, ...
          -0.060, -0.058, -0.062, -0.059;
           0.009, -0.011, -0.023, -0.037, -0.050, -0.047, -0.074, -0.079, ...
          -0.091, -0.076, -0.077, -0.076];

% 2-D Interpolation.
cli = interp2(alphtab_deg, ...
              betatab_deg, ...
              cltabl, ...
              alpha_ltd_deg, ...
              abs(beta_ltd_deg));

% Account for the sign of the sideslip.
if beta_deg < 0.0
    cli = -cli;
end