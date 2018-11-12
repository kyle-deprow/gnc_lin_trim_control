function tgear = f16_tgear(thtl)
% Data and interpolation for F-16 power command vs. throttle.
% 
% Inputs:
%   thtl   scalar, throttle position (decimal percentage)
%
% Outputs:
%   tgear  scalar, power command
%
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

if thtl <= 0.77
    tgear = 64.94*thtl;
else
    tgear = 217.38*thtl - 117.38;
end
