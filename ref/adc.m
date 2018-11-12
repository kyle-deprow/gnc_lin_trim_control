function [mach, qbar_psf] = adc(v_fps, h_ft)
% This function estimates the mach number and dynamic pressure for a
% particular flight condition.
%
%
% Inputs:
%   v_fps       scalar, velocity (fps)
%   h_ft        scalar, altitude (ft)
%
%
% Outputs:
%   mach        scalar, mach number
%   qbar_psf    scalar, dynamic pressure (psf)
%
% 
% Version history:
%   September 6, 2011

% Parameters
% Sea level density, R0.
R0_cf = 2.377e-3;

% Estimate the temperature.
tfac = 1.0 - 0.703e-5 * h_ft;
if h_ft >= 35000.0
    T = 390;
else
    T = 519.0*tfac;
end

% Estimate the air density at altitude.
rho_cf = R0_cf * (tfac^4.14);

% Estimate Mach number.
mach = v_fps/sqrt(1.4*1716.3*T);

% Estimate dynamic pressure.
qbar_psf = 0.5 * rho_cf * v_fps^2;
