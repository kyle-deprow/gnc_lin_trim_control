function pdot = f16_del_power(p3, p1)
% F-16 power rate of change
% 
%
% Inputs:
%   p3      scalar, actual power
%   p1      scalar, power command
%
%
% Outputs:
%   pdot    scalar, power change
%
%
% September 8, 2011
%
%
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

if p1 >= 50
    if p3 >= 50
        T = 5.0;
        p2 = p1;
    else
        p2 = 60;
        T = rtau(p2-p3);
    end
else
    if p3 >= 50
        T = 5.0;
        p2 = 40.0;
    else
        p2 = p1;
        T = rtau(p2-p3);
    end
end

pdot = T*(p2-p3);


function T = rtau(dp)
% Computes the reciprocal time constant.
if dp <= 25
    T = 1.0;
elseif dp >= 50
    T = 0.1;
else
    T = 1.9 - 0.036*dp;
end
