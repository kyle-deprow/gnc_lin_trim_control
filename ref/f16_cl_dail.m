function clda = f16_cl_dail(alpha_deg, beta_deg)
% Data and interpolation for F-16 Rolling Moment due to Aileron deflection
% Coefficient.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip, in degrees
%
% Outputs:
%   clda     rolling moment due to aileron deflection coefficient
% 
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Independent variables' breakpoints.
alphtab_deg = -10:5:45;
betatab_deg = -15:5:15;

% Bound alpha and beta to the table range.
alpha_ltd_deg = max(-10, min(45, alpha_deg));
beta_ltd_deg  = max(-15, min(15, beta_deg));

% Clda table
clda_tab=[-0.041, -0.052, -0.053, -0.056, -0.050, -0.056, -0.082, -0.059, ...
          -0.042, -0.038, -0.027, -0.017;
          -0.041, -0.053, -0.053, -0.053, -0.050, -0.051, -0.066, -0.043, ...
          -0.038, -0.027, -0.023, -0.016;
          -0.042, -0.053, -0.052, -0.051, -0.049, -0.049, -0.043, -0.035, ...
          -0.026, -0.016, -0.018, -0.014;
          -0.040, -0.052, -0.051, -0.052, -0.048, -0.048, -0.042, -0.037, ...
          -0.031, -0.026, -0.017, -0.012;
          -0.043, -0.049, -0.048, -0.049, -0.043, -0.042, -0.042, -0.036, ...
          -0.025, -0.021, -0.016, -0.011;
          -0.044, -0.048, -0.048, -0.047, -0.042, -0.041, -0.020, -0.028, ...
          -0.013, -0.014, -0.011, -0.010;
          -0.043, -0.049, -0.047, -0.045, -0.042, -0.037, -0.003, -0.013, ...
          -0.010, -0.003, -0.007, -0.008];

% 2-D Interpolation.
clda = interp2(alphtab_deg, ...
               betatab_deg, ...
               clda_tab, ...
               alpha_ltd_deg, ...
               beta_ltd_deg, ...
               '*linear');
