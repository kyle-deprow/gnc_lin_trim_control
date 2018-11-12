function cxi = f16_cx(alpha_deg, elev_deg)
% Data and interpolation for F-16 Axial Force Coefficient.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   elev_deg    elevator deflection, in degrees
%
% Outputs:
%   cxi     axial force coefficient
% 
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Independent variables' breakpoints.
alphtab_deg = -10:5:45;
elevtab_deg = -24:12:24;

% Bound aoa and elevator to table ranges.
alpha_ltd_deg = max(-10, min(45, alpha_deg));
elev_ltd_deg  = max(-24, min(24, elev_deg));

% CX table
cxtabl = [-0.099, -0.081, -0.081, -0.063, -0.025, 0.044, 0.097, 0.113, ...
          0.145, 0.167, 0.174, 0.166;
          -0.048, -0.038, -0.040, -0.021, 0.016, 0.083, 0.127, 0.137, ...
          0.162, 0.177, 0.179, 0.167;
          -0.022, -0.020, -0.021, -0.004, 0.032, 0.094, 0.128, 0.130, ...
          0.154, 0.161, 0.155, 0.138;
          -0.040, -0.038, -0.039, -0.025, 0.006, 0.062, 0.087, 0.085, ...
          0.100, 0.110, 0.104, 0.091;
          -0.083, -0.073, -0.076, -0.072, -0.046, 0.012, 0.024, 0.025, ...
          0.043, 0.053, 0.047, 0.040];
      
% 2-D Interpolation.
cxi = interp2(alphtab_deg, ...
              elevtab_deg, ...
              cxtabl, ...
              alpha_ltd_deg, ...
              elev_ltd_deg, ...
              '*linear');
