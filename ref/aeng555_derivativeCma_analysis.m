% AENG-555: Guidance and Control of Aerospace Vehicles
% Saint Louis University
% Ken Buckholtz, Ph.D.
%
% October 25, 2018
% This script determines the derivative of Cmalpha using discrete values of alphas


clear 
close all

% Define constants.
d2r = pi/180;

% Define range of AOA to plot over.
alpha_deg = 0:.1:30;
n_alpha   = numel(alpha_deg);

% Define Elevator breakpoints to plot over.
dele_deg = [0];
n_dele   = numel(dele_deg);

% Define sideslip angle.
beta_deg = 0;

% Pre-allocate memory.
% Create vectors to store the various coefficients.
Cmalpha = zeros(n_dele,n_alpha);

% Define sideslip angle.
beta_deg = 0;

% Loop over each elevator channel and compute coefficients at each angle of
% attack.
for i_dele = 1:n_dele
    for i_alpha = 1:n_alpha
        % Get the pitching moment coefficient.
        Cmalpha(i_dele,i_alpha) = f16_cm(alpha_deg(i_alpha), ...
                                    dele_deg(i_dele));
    end
end

% Generate first order derivative
% find h, the derivative interval
h = (alpha_deg(end) - alpha_deg(1))/size(alpha_deg,2);
% Call first_order_derivative
derivative_Cmalpha_alpha = first_order_derivative(Cmalpha, h);

figure(1)
plot(alpha_deg, Cmalpha); grid on;
xlabel('\alpha deg');
ylabel('C_m\alpha');
legend('\delta_e = 0 deg');
title('C_m\alpha vs \alpha');

figure(2)
plot(alpha_deg, derivative_Cmalpha_alpha); grid on;
xlabel('\alpha deg');
ylabel('derivative dC_m_\alpha/\alpha');
legend('\delta_e = 0 deg');
title('Derivative dC_m_\alpha/d\alpha vs \alpha');
