% AENG-555: Guidance and Control of Aerospace Vehicles
% Saint Louis University
% Kyle DeProw
%
% October 25, 2018 
% This script plots thrust over an alpha range and over varius
% elevator actuations


clear 
close all

% Define power range.
power = 0:10:100;
n_power   = numel(power);

% Define altitude range
h_ft = [25000, 50000];
n_h   = numel(h_ft);

% define Mach range
mach = [0.4:0.1:2];
n_mach = numel(mach);

% Define sideslip angle.
beta_deg = 0;

% Pre-allocate memory.
% Create vectors to store the various coefficients.
thrust_power = zeros(n_h,n_power);
thrust_mach = zeros(n_h,n_mach);

for i_h = 1:n_h
    for i_power = 1:n_power
        % Get the thrust with varying power and h at mach=06
        thrust_power(i_h,i_power) = f16_thrust(power(i_power), ...
                                        h_ft(i_h), 0.6);
    end
end

for i_h = 1:n_h
    for i_mach = 1:n_mach
        % Get the thrust with varying mach and h at power = 100%
        thrust_mach(i_h,i_mach) = f16_thrust(100, ...
                                        h_ft(i_h), mach(i_mach));
    end
end

figure(1)
plot(power, thrust_power); grid on;
xlabel('power %');
ylabel('thrust');
legend('h = 25kft', 'h = 50kft');
title('Thrust vs power at various altitudes');

figure(2)
plot(mach, thrust_mach); grid on;
xlabel('Mach Number');
ylabel('thrust (lbs)');
legend('h = 25kft', 'h = 50kft');
title('Thrust vs Mach at various altitudes');
