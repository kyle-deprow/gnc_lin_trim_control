function cldr = f16_cl_drud(alpha_deg, beta_deg)
% Data and interpolation for F-16 Rolling Moment due to Rudder deflection
% Coefficient.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip, in degrees
%
% Outputs:
%   cldr     rolling moment due to rudder deflection coefficient
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
cldr_tab=[ 0.005,  0.017,  0.014,  0.010, -0.005,  0.009,  0.019,  0.005, ...
          -0.000, -0.005, -0.011,  0.008;
           0.007,  0.016,  0.014,  0.014,  0.013,  0.009,  0.012,  0.005, ...
           0.000,  0.004,  0.009,  0.007;
           0.013,  0.013,  0.011,  0.012,  0.011,  0.009,  0.008,  0.005, ...
          -0.002,  0.005,  0.003,  0.005;
           0.018,  0.015,  0.015,  0.014,  0.014,  0.014,  0.014,  0.015, ...
           0.013,  0.011,  0.006,  0.001;
           0.015,  0.014,  0.013,  0.013,  0.012,  0.011,  0.011,  0.010, ...
           0.008,  0.008,  0.007,  0.003;
           0.021,  0.011,  0.010,  0.011,  0.010,  0.009,  0.008,  0.010, ...
           0.006,  0.005,  0.000,  0.001;
           0.023,  0.010,  0.011,  0.011,  0.011,  0.010,  0.008,  0.010, ...
           0.006,  0.014,  0.020,  0.000];

% 2-D Interpolation.
cldr = interp2(alphtab_deg, ...
               betatab_deg, ...
               cldr_tab, ...
               alpha_ltd_deg, ...
               beta_ltd_deg, ...
               '*linear');
