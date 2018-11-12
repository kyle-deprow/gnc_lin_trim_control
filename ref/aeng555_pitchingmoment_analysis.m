% AENG-555: Guidance and Control of Aerospace Vehicles
% Saint Louis University
% Kyle DeProw
%
% October 25, 2018
% This script plots pitching moment coefficient over an alpha
% range and over varius elevator actuations


clear 
close all

% Define constants.
d2r = pi/180;

% Define range of AOA to plot over.
alpha_deg = -10:0.1:50;
n_alpha   = numel(alpha_deg);

% Define Elevator breakpoints to plot over.
dele_deg = [-25, -10, 0, 10, 25];
n_dele   = numel(dele_deg);

% Define sideslip angle.
beta_deg = 0;

% Pre-allocate memory.
% Create vectors to store the various coefficients.
CM = zeros(n_dele,n_alpha);

% Loop over each elevator channel and compute coefficients at each angle of
% attack.
for i_dele = 1:n_dele
    for i_alpha = 1:n_alpha
        % Get the pitching moment coefficient.
        CM(i_dele,i_alpha) = f16_cm(alpha_deg(i_alpha), ...
                                    dele_deg(i_dele));
    end
end

figure(1)
plot(alpha_deg, CM); grid on;
xlabel('\alpha, deg');
ylabel('CM');
legend('\delta_e = -25 deg', '\delta_e = -10 deg', '\delta_e = 0 deg','\delta_e = 10 deg', '\delta_e = 25 deg');
title('Pitching moment coefficient, Cm, vs. \alpha for various elevator deflections');

