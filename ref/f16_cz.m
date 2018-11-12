function czi = f16_cz(alpha_deg, beta_deg, elev_deg)
% Data and interpolation for F-16 z-Force Coefficient.
% 
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip angle, in degrees
%   elev_deg    elevator deflection, in degrees
%
% Outputs:
%   cxi     axial force coefficient
%
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Constant.
d2r = pi/180;

% Scale the input to work with data table.
%alpha_deg = alpha_deg/5 + 3;

% Limit the elevator deflection to 25 deg.
elev_ltd_deg = min(25,max(-25,elev_deg));

alphtab_deg = -10:5:45;

cztabl = [0.770, 0.241, -0.100, -0.416, -0.731, -1.053, -1.366, -1.646, ...
          -1.917, -2.120, -2.248, -2.229];
      
% 1-D Interpolation.
c1 = interp1(alphtab_deg, cztabl, alpha_deg);

% z-Force calculation.
czi = c1 * (1-(beta_deg*d2r)^2) - 0.19*(elev_ltd_deg/25.0);