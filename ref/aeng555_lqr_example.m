%**************************************************************************
%
% AENG-555 : Automatic Control Systems
% Saint Louis University
% Ken Buckholtz, Ph.D.
%
% Continuing MATLAB example following the flow as laid out in the text:
%   Linear State-Space Control Systems
%   by Williams II and Lawrence, Wiley, 2007
%
% The system under study is a complex mechanical system involving two
% rotating coupled inertias and a vertical translating mass.  The input
% into the system is a torque applied to a body.  The outputs are the
% angular displacments of the inertias and the tranlational displacement of
% the mass.
%
% Version history:
%   2016.11.02  KRB  Adapted from AENG 556 Class
%
%**************************************************************************

set(0,'DefaultTextInterpreter','tex');
clear
close all

dtr = pi/180;
rtd = 180/pi;

% Step command used for step response analysis.
step_mag = 1*dtr;

%%%%%%%%%% State Space Model %%%%%%%%%%
% Define the state space model for the three body problem.
%sys_ss = ContExample_model();
% Parameter definitions
%M   = 10;   % kg
%B1  = 10;   % Nm-sec
%J1  = 15;   % Nm-sec^2
%J2  = 15;   % Nm-sec^2
%B2  = 100;  % N-sec/m
%K1  = 100;  % Nm/rad
%K2  = 130;  % N/m
%K3  = 160;  % N/m
%R   = 1;    % m

%% Define the state space matrices (A, B).
%A   = [0, 1, 0, 0, 0, 0;
       %-K1/J1, 0, K1/J1, 0, 0, 0;
       %0, 0, 0, 1, 0, 0;
       %K1/J2, 0, -(K1+R^2*K2)/J2, -B1/J2, R*K2/J2, 0;
       %0, 0, 0, 0, 0, 1;
       %0, 0, R*K2/M, 0, -(K2+K3)/M, -B2/M];
   
%B   = [0; 1/J1; 0; 0; 0; 0];

A = [-0.008467, 37.4, -32.17, 0.7555; -0.00049587, -0.485739, 0, 0.9244; 0,0,0,1; 0, -0.1832, 0, -0.41797];
B = [0.0547; -0.00091; 0; -0.03668];

[n,m] = size(B);

% Define the output matrices (C, D).
C   = eye(n);

D   = zeros(n,m);

% Define the state space model.
sys_ss = ss(A,B,C,D);

% Define system size parameters.
p = size(sys_ss.c,1);


%%%%%%%%%% Modal Analysis %%%%%%%%%%
% Compute the sensitivities and metrics of the modes of the system.
[Msens, Mmetrics] = perform_mode_analysis(A);

    
%s%%%%%%%%% Optimal Control: Tracker %%%%%%%%%%
% Define to track rotational position of J2.
Creg = [0,0,1,0];
Dreg = sys_ss.d(3,:);
nreg = 1;

% Build the servomechanizm state space model of the plant by appending an
% integral error state for each regulating channel.
Aw = zeros(n+nreg, n+nreg);
Aw(1:nreg,nreg+1:end) = Creg;
Aw(nreg+1:end,nreg+1:end) = sys_ss.a;

Bw = zeros(n+nreg, m);
Bw(1:nreg,:) = Dreg;
Bw(nreg+1:end,:) = sys_ss.b;


% Specify range of scale factor to apply to the integral error state.
q = logspace(1,6,300);

% Specify the LQR matrices.
Qw = zeros(n+nreg, n+nreg);
Qw(n+nreg,n+nreg) = 3000;

Rw = 1*eye(m);
Rw = 0.01*eye(m);

% Define frequency range for performing frequency analysis.
w = logspace(-3,3,600);

% Form plant system.  Set output matrix for state feedback control.
sys_p.Ap = sys_ss.A;
sys_p.Bp = sys_ss.B;
sys_p.Cp = eye(n);
sys_p.Dp = zeros(n,1);

% Pass in the characteristics of the system.  Note, this process % returns plots.  This process assumes Dreg == 0.
control_data = perform_charting(sys_p, Aw, Bw, Creg, Qw, Rw, q, w);
    

%%%%%%%%%% Perform Selected Design %%%%%%%%%%
% Set the integral error state weights to the desired scale factor.
for j = 1:nreg
    Qw(j,j) = 400000;
end

% Solve the LQR problem.
[Kw, S, E] = lqr(Aw, Bw, Qw, Rw);
Kw

% Set the state feedback gain and the integral error gain.
kI = Kw(:,1:nreg);
Kv = Kw(:,nreg+1:end);

% Build controller object to track reference theta2.
sys_c.Ac = zeros(nreg,nreg);
sys_c.Bc1= Creg;
sys_c.Bc2= -eye(nreg);
sys_c.Cc = -kI;
sys_c.Dc1= -Kv;
sys_c.Dc2= zeros(m,nreg);

% Form closed loop system.
[Acl, Bcl, Ccl, Dcl] = closed_loop_system(sys_p, sys_c);

%% Perform a step response.
syscl = ss(Acl,Bcl,Ccl,Dcl);
[yoptimal, timeoptimal_sec, x] = step(syscl*0.1,3);

%% Reconstruct the control.
uoptimal = -Kw(:,[nreg+1:nreg+n,1:nreg])*x';

figure
plot(timeoptimal_sec,yoptimal(:,3)*rtd,'-b', 'LineWidth', 2);grid on
xlabel('Time, sec');
ylabel('\theta, deg');
title('Step response for servomechanism designs');


figure
plot(timeoptimal_sec, uoptimal,'-b', 'LineWidth', 2); grid on;
xlabel('Time, sec');
ylabel('\delta_{Elevator}, degrees');
title('Control for servomechanism designs');


%%%%% FREQUENCY ANALYSIS %%%%%
% Perform frequency analysis of servo tracker.
%freq_resp = freq_analysis(sys_p, sys_c, w, 1, 3);

%% Plot frequency data.
%plot_freq_resp(freq_resp, w);
