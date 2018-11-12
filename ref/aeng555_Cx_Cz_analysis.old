% AENG-555: Guidance and Control of Aerospace Vehicles
% Saint Louis University
% Ken Buckholtz, Ph.D.
%
% July 22, 2016
% This script plots logitudinal forces in body axis system across elevator
% settings.


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
CX = zeros(n_dele,n_alpha);
CZ = CX;

% Loop over each elevator channel and compute coefficients at each angle of
% attack.
for i_dele = 1:n_dele
    for i_alpha = 1:n_alpha
        % Get the axial coefficient.
        CX(i_dele,i_alpha) = f16_cx(alpha_deg(i_alpha), ...
                                    dele_deg(i_dele));

        % Get the vertical coefficient.
        CZ(i_dele,i_alpha) = f16_cz(alpha_deg(i_alpha), ...
                                    beta_deg, ...
                                    dele_deg(i_dele));
    
        % Calculate Lift coefficient
        CL(i_dele,i_alpha) = sin(alpha_deg(i_alpha)*d2r)*CX(i_dele,i_alpha) - cos(alpha_deg(i_alpha)*d2r)*CZ(i_dele,i_alpha);
    end
end



% Generate plots.
figure(1)
plot(alpha_deg, CL); grid on;
xlabel('\alpha, deg');
ylabel('Cx');
legend('\delta_e = -25 deg', '\delta_e = -10 deg', '\delta_e = 0 deg','\delta_e = 10 deg', '\delta_e = 25 deg');
title('Lift coefficient, CL, vs. \alpha for various elevator deflections');

