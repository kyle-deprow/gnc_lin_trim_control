function cndr = f16_cn_drud(alpha_deg, beta_deg)
% Data and interpolation for F-16 Yawing Moment due to Rudder deflection
% Coefficient.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip, in degrees
%
% Outputs:
%   cndr     yawing moment due to rudder deflection coefficient
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

% Cndr table
cndr_tab=[-0.018, -0.052, -0.052, -0.052, -0.054, -0.049, -0.059, -0.051, ...
          -0.030, -0.037, -0.026, -0.013;
          -0.028, -0.051, -0.043, -0.046, -0.045, -0.049, -0.057, -0.052, ...
          -0.030, -0.033, -0.030, -0.008;
          -0.037, -0.041, -0.038, -0.040, -0.040, -0.038, -0.037, -0.030, ...
          -0.027, -0.024, -0.019, -0.013;
          -0.048, -0.045, -0.045, -0.045, -0.044, -0.045, -0.047, -0.048, ...
          -0.049, -0.045, -0.033, -0.016;
          -0.043, -0.044, -0.041, -0.041, -0.040, -0.038, -0.034, -0.035, ...
          -0.035, -0.029, -0.022, -0.009;
          -0.052, -0.034, -0.036, -0.036, -0.035, -0.028, -0.024, -0.023, ...
          -0.020, -0.016, -0.010, -0.014;
          -0.062, -0.034, -0.027, -0.028, -0.027, -0.027, -0.023, -0.023, ...
          -0.019, -0.009, -0.025, -0.010];

% 2-D Interpolation.
cndr = interp2(alphtab_deg, ...
               betatab_deg, ...
               cndr_tab, ...
               alpha_ltd_deg, ...
               beta_ltd_deg, ...
               '*linear');
