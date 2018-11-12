function damp_coeffs = f16_damp(alpha_deg)
% Function to compute the various damping derivatives for the F-16.
%
% Inputs:
%   alpha_deg   angle of attack, in degrees
%   beta_deg    sideslip, in degrees
%
% Outputs:
%   damp_coeffs     structure, contains the computed damping coefficients.
%       .CXq
%       .CYr
%       .CYp
%       .CZq
%       .Clr
%       .Clp
%       .Cmq
%       .Cnr
%       .Cnp
% 
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Independent variables' breakpoints.
alphtab_deg = -10:5:45;

% Damping table.
damp_tab =[ -0.267, -0.110,  0.308,  1.340,  2.080,  2.910,  2.760,  2.050, ...
             1.500,  1.490,  1.830,  1.210;
             0.882,  0.852,  0.876,  0.958,  0.962,  0.974,  0.819,  0.483, ...
             0.590,  1.210, -0.493, -1.040;
            -0.108, -0.108, -1.880,  0.110,  0.258,  0.226,  0.344,  0.362, ...
             0.611,  0.529,  0.298, -2.270;
            -8.800,-25.800,-28.900,-31.400,-31.200,-30.700,-27.700,-28.200,...
           -29.000,-29.800,-38.300,-35.300;
            -0.126, -0.026,  0.063,  0.113,  0.208,  0.230,  0.319,  0.437,...
             0.680,  0.100,  0.447, -0.330;
            -0.360, -0.359, -0.443, -0.420, -0.383, -0.375, -0.329, -0.294,...
            -0.230, -0.210, -0.120, -0.100;
            -7.210, -0.540, -5.230, -5.260, -6.110, -6.640, -5.690, -6.000,...
            -6.200, -6.400, -6.600, -6.000;
            -0.380, -0.363, -0.378, -0.386, -0.370, -0.453, -0.550, -0.582, ...
            -0.595, -0.637, -1.020, -0.840;
             0.061,  0.052,  0.052, -0.012, -0.013, -0.024,  0.050,  0.150, ...
             0.130,  0.158,  0.240,  0.150];
         
% Loop through computing the damping coefficients.
d_coeff = zeros(1,9);
for i = 1:9
    d_coeff(i) = interp1(alphtab_deg, ... 
                         damp_tab(i,:), ...
                         alpha_deg, ...
                         '*linear');
end
           

% Store the coefficients.
damp_coeffs.CXq = d_coeff(1);
damp_coeffs.CYr = d_coeff(2);
damp_coeffs.CYp = d_coeff(3);
damp_coeffs.CZq = d_coeff(4);
damp_coeffs.Clr = d_coeff(5);
damp_coeffs.Clp = d_coeff(6);
damp_coeffs.Cmq = d_coeff(7);
damp_coeffs.Cnr = d_coeff(8);
damp_coeffs.Cnp = d_coeff(9);