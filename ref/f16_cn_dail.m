function cnda = f16_cn_dail(alpha_deg, beta_deg)
% Data and interpolation for F-16 Yawing Moment due to Aileron deflection
% Coefficient.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip, in degrees
%
% Outputs:
%   cnda     yawing moment due to aileron deflection coefficient
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

% Cnda table
cnda_tab=[ 0.001, -0.027, -0.017, -0.013, -0.012, -0.016,  0.001,  0.017, ...
           0.011,  0.017,  0.008,  0.016;
           0.002, -0.014, -0.016, -0.016, -0.014, -0.019, -0.021,  0.002, ...
           0.012,  0.015,  0.015,  0.011;
          -0.006, -0.008, -0.006, -0.006, -0.005, -0.008, -0.005,  0.007, ...
           0.004,  0.007,  0.006,  0.006;
          -0.011, -0.011, -0.010, -0.009, -0.008, -0.006,  0.000,  0.004, ...
           0.007,  0.010,  0.004,  0.010;
          -0.015, -0.015, -0.014, -0.012, -0.011, -0.008, -0.002,  0.002, ...
           0.006,  0.012,  0.011,  0.011;
          -0.024, -0.010, -0.004, -0.002, -0.001,  0.003,  0.014,  0.006, ...
          -0.001,  0.004,  0.004,  0.006;
          -0.022,  0.002, -0.003, -0.005, -0.003, -0.001, -0.009, -0.009, ...
          -0.001,  0.003, -0.002,  0.001];

% 2-D Interpolation.
cnda = interp2(alphtab_deg, ...
               betatab_deg, ...
               cnda_tab, ...
               alpha_ltd_deg, ...
               beta_ltd_deg, ...
               '*linear');
