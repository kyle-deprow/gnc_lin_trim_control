function ti_lbs = f16_thrust(pow, alt_ft, mach)
% Data and interpolation for F-16 thrust.
% 
% Inputs:
%   pow     percentage of power
%   alt_ft  aircraft altitude (in feet)
%   mach    Mach number
%
% Outputs:
%   ti_lbs  thrust (in lbs)
%
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Independent variable breakpoints.
htabl_ft = 0:1e4:5e4;
machtabl = 0:0.2:1;

% Thrust tables depending upon power.
idle = [1060, 670, 880, 1140, 1500, 1860;
        635, 425, 690, 1010, 1330, 1700;
        60, 25, 345, 755, 1130, 1525;
        -1020, -710, -300, 350, 910, 1360;
        -2700, -1900, -1300, -247, 600, 1100;
        -3600, -1400, -595, -342, -200, 700];
mil = [12680, 9150, 6200, 3950, 2450, 1400;
       12680, 9150, 6313, 4040, 2470, 1400;
       12610, 9312, 6610, 4290, 2600, 1560;
       12640, 9839, 7090, 4660, 2840, 1660;
       12390, 10176, 7750, 5320, 3250, 1930;
       11680, 9848, 8050, 6100, 3800, 2310];
maxtab =   [20000, 15000, 10800, 7000, 4000, 2500;
            21420, 15700, 11225, 7323, 4435, 2600;
            22700, 16860, 12250, 8154, 5000, 2835;
            24240, 18910, 13760, 9285, 5700, 3215;
            26070, 21075, 15975, 11115, 6860, 3950;
            28886, 23319, 18300, 13484, 8642, 5057];
% Get thrust from the military table.
milthr = interp2(htabl_ft, machtabl, mil, alt_ft, mach);
if pow < 50
    % Compute thrust between idle and military tables.
    idlethr = interp2(htabl_ft, machtabl, idle, alt_ft, mach);
    ti_lbs  = idlethr + (milthr - idlethr)*pow*0.02;
else
    % Compute thrust between military and max tables.
    maxthr = interp2(htabl_ft, machtabl, maxtab, alt_ft, mach);
    ti_lbs = milthr + (maxthr - milthr)*(pow-50)*0.02;
end
