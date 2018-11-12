function [xd, Az, Ay] = f16_nonlinear_model(time, x, u)
% Nonlinear 6DOF aircraft model.
%
% Inputs:
%   time    scalar, instance of time (dummy)
%   x       n element vector, state vector
%       x(1)... True Airspeed, Vt (ft/sec)
%       x(2)... Angle of Attack, alpha (rad)
%       x(3)... Angle of Sideslip, beta (rad)
%       x(4)... Roll attitude, phi (rad)
%       x(5)... Pitch attitude, theta (rad)
%       x(6)... Yaw Attitude, psi (rad)
%       x(7)... Roll rate, P (rad/sec)
%       x(8)... Pitch rate, Q (rad/sec)
%       x(9)... Yaw rate, R (rad/sec)
%       x(10).. North position, N (ft)
%       x(11).. East position, E (ft)
%       x(12).. Altitude, h (ft)
%       x(13).. power (W)
%   u       m element vector, control vector
%       u(1)... throttle, thtl (0-1)
%       u(2)... Elevator, de (deg)
%       u(3)... Aileron, da (deg)
%       u(4)... Rudder, dr (deg)
%
%
% Outputs:
%   xd      n element vector, state derivative vector
%   Az      scalar, vertical acceleartion (g)
%   Ay      scalar, lateral acceleration (g)
%
% 
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis

% Define globals.
global xcg

% Parameters
axx = 9496.0;
ayy = 55814.0;
azz = 63100.0;
axz = 982.0;
weight_lbs = 20490.466; % 25000.0;
s_ft2 = 300;
b_ft = 30;
cbar_ft = 11.32;
xcgr = 0.35;
hx_slgft2 = 160.0;

% Initialize the output
xd = zeros(size(x));


% Initialize additional parameters and constants.
g_fps2  = 32.174;
mass_slg = weight_lbs/g_fps2;
r2d     = 180/pi;

axzs    = axz^2;
xpq     = axz*(axx-ayy+azz);
gam     = axx*azz - axz^2;
xqr     = azz*(azz-ayy)+axzs;
zpq     = (axx-ayy)*axx+axzs;
ypr     = azz - axx;


% Assign states to local variables.  Convert from radians to degress for
% AOA and AOS.
vt_fps      = x(1);
alpha_deg   = x(2)*r2d;
beta_deg    = x(3)*r2d;
phi_rad     = x(4);
theta_rad   = x(5);
psi_rad     = x(6);
P_rps       = x(7);
Q_rps       = x(8);
R_rps       = x(9);
alt_ft      = x(12);
pow         = x(13);


% Assign controls to local variables.  Leave surfaces in terms of degrees.
thtl    = u(1);
de_deg  = u(2);
da_deg  = u(3);
dr_deg  = u(4);


% Compute air data parameters.
[mach,...
 qbar_psf] = adc(vt_fps, ...
                 alt_ft);
 
             
             
% Extract data from the engine model.
cpow = f16_tgear(thtl);
xd(13) = f16_del_power(pow,cpow);
T_lbs = f16_thrust(pow, alt_ft, mach);


% Aerodynamic coefficient lookups.
cxt = f16_cx(alpha_deg, de_deg);
cyt = f16_cy(beta_deg, da_deg, dr_deg);
czt = f16_cz(alpha_deg, beta_deg, de_deg);

clt = f16_cl(alpha_deg, beta_deg) + ...
      f16_cl_dail(alpha_deg, beta_deg)*(da_deg/20) + ...
      f16_cl_drud(alpha_deg, beta_deg)*(dr_deg/30);
cmt = f16_cm(alpha_deg, de_deg);

cnt = f16_cn(alpha_deg, beta_deg) + ...
      f16_cn_dail(alpha_deg, beta_deg)*(da_deg/20) + ...
      f16_cn_drud(alpha_deg, beta_deg)*(dr_deg/30);
  

% Damping derivatives.
tvt_s = 0.5/vt_fps;
b2v = b_ft*tvt_s;
cq = cbar_ft*Q_rps*tvt_s;

damp_deriv = f16_damp(alpha_deg);
cxt = cxt + cq * damp_deriv.CXq;
cyt = cyt + b2v * (damp_deriv.CYr*R_rps + damp_deriv.CYp*P_rps);
czt = czt + cq * damp_deriv.CZq;
clt = clt + b2v * (damp_deriv.Clr*R_rps + damp_deriv.Clp*P_rps);
cmt = cmt + cq * damp_deriv.Cmq + czt * (xcgr - xcg);
cnt = cnt + b2v * (damp_deriv.Cnr*R_rps + damp_deriv.Cnp*P_rps) - ...
      cyt * (xcgr - xcg) * (cbar_ft / b_ft);
  
  
% Pre-compute some variables for the state space model.  Computations for
% AOA and AOS must be in radians.
salp = sin(x(2));
calp = cos(x(2));
sbta = sin(x(3));
cbta = cos(x(3));
sth = sin(theta_rad);
cth = cos(theta_rad);
sph = sin(phi_rad);
cph = cos(phi_rad);
spsi= sin(psi_rad);
cpsi= cos(psi_rad);

qs = qbar_psf * s_ft2;
qsb = qs * b_ft;
rmqs = qs/mass_slg;
gcth = g_fps2 * cth;
qsph = Q_rps * sph;
Ay_fps2 = rmqs * cyt;
Az_fps2 = rmqs * czt;

% Velocities along wind axes.
U_fps = vt_fps * calp*cbta;
V_fps = vt_fps * sbta;
W_fps = vt_fps * salp*cbta;


% Force equations.
udot_fps2 = R_rps * V_fps - Q_rps * W_fps - g_fps2 * sth + ...
            (qs * cxt + T_lbs) / mass_slg;
vdot_fps2 = P_rps * W_fps - R_rps * U_fps + gcth * sph + Ay_fps2;
wdot_fps2 = Q_rps * U_fps - P_rps * V_fps + gcth * cph + Az_fps2;
dum = (U_fps^2 + W_fps^2);

xd(1) = (U_fps*udot_fps2 + V_fps*vdot_fps2 + W_fps*wdot_fps2) / vt_fps;
xd(2) = (U_fps*wdot_fps2 - W_fps*udot_fps2) / dum;
xd(3) = (vt_fps*vdot_fps2 - V_fps*xd(1)) *cbta / dum;

% Kinematics.
xd(4) = P_rps + (sth/cth) * (qsph + R_rps * cph);
xd(5) = Q_rps * cph - R_rps * sph;
xd(6) = (qsph + R_rps * cph) / cth;

% Moments.
roll = qsb * clt;
pitch = qs * cbar_ft * cmt;
yaw = qsb * cnt;
PQ = P_rps * Q_rps;
QR = Q_rps * R_rps;
PR = P_rps * R_rps;
QHX = Q_rps * hx_slgft2;

xd(7) = ( xpq*PQ - xqr*QR + azz*roll + axz*(yaw + QHX) ) / gam;
xd(8) = ( ypr*PR - axz*(P_rps^2 - R_rps^2) + pitch - R_rps*hx_slgft2) / ayy;
xd(9) = ( zpq*PQ - xpq*QR + axz*roll + axx*(yaw + QHX)) / gam;

% Navigation.
t1 = sph*cpsi;
t2 = cph*sth;
t3 = sph*spsi;
s1 = cth*cpsi;
s2 = cth*spsi;
s3 = t1*sth - cph*spsi;
s4 = t3*sth + cph*cpsi;
s5 = sph*cth;
s6 = t2*cpsi + t3;
s7 = t2*spsi - t1;
s8 = cph*cth;

xd(10) = U_fps*s1 + V_fps*s3 + W_fps*s6;
xd(11) = U_fps*s2 + V_fps*s4 + W_fps*s7;
xd(12) = U_fps*sth - V_fps*s5 - W_fps*s8;

% Set the outputs of the vertical and lateral accelerations, in g's
Az = Az_fps2/g_fps2;
Ay = Ay_fps2/g_fps2;
