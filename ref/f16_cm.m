function cmi = f16_cm(alpha_deg, elev_deg)
% Data and interpolation for F-16 Pitching Moment.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   elev_deg    elevator deflection, in degrees
%
% Outputs:
%   cmi    pitching moment coefficient
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

% CM table.
cmtabl = [ 0.205,  0.168,  0.186,  0.196,  0.213,  0.251,  0.245,  0.238, ...
           0.252,  0.231,  0.198,  0.192;
           0.081,  0.077,  0.107,  0.110,  0.110,  0.141,  0.127,  0.119, ...
           0.133,  0.108,  0.081,  0.093;
          -0.046, -0.020, -0.009, -0.005, -0.006,  0.010,  0.006, -0.001, ...
           0.014,  0.000, -0.013,  0.032;
          -0.174, -0.145, -0.121, -0.127, -0.129, -0.102, -0.097, -0.113, ...
          -0.087, -0.084, -0.069, -0.006;
          -0.259, -0.202, -0.184, -0.193, -0.199, -0.150, -0.160, -0.167,...
          -0.104, -0.076, -0.041, -0.005];
      
 % 2-D interpolation.
 cmi = interp2(alphtab_deg, ...
               elevtab_deg, ...
               cmtabl, ...
               alpha_ltd_deg, ...
               elev_ltd_deg, ...
               '*linear');