function [A,B,C,D] = closed_loop_system(sys_p, sys_c)
%**************************************************************************
%
% Function computes the closed loop matrices for a particular plant and
% controller.
%
%
% Inputs:
%   sys_p   structure, plant state space system
%           Fields:
%               .Ap     [nxn] matrix, State matrix
%               .Bp     [nxm] matrix, Control matrix
%               .Cp     [pxn] matrix, Output matrix
%               .Dp     [pxm] matrix, Direct transmission matrix
%   sys_c   structure, controller state space system
%           Fields:
%               .Ac     [ncxnc] matrix, Controller state matrix
%               .Bc1    [ncxmc1] matrix, Plant output matrix
%               .Bc2    [ncxnr] matrix, Reference matrix
%               .Cc     [mxnc] matrix, Controller output matrix
%               .Dc1    [mxmc1] matrix, Control-Plant output matrix
%               .Dc2    [mxnr] matrix, Control-Reference matrix
%
%
% Outputs:
%   A   [nxn] matrix, Closed loop state matrix
%   B   [nxm] matrix, Closed loop control matrix
%   C   [pxn] matrix, Closed loop output matrix
%   D   [pxm] matrix, Closed loop direct transmission matrix
%
%
% Version history:
%   2013.12.30  KRB  Initial release
%
%**************************************************************************

% Extract plant matrices.
Ap = sys_p.Ap;
Bp = sys_p.Bp;
Cp = sys_p.Cp;
Dp = sys_p.Dp;

% Extract controller matrices.
Ac = sys_c.Ac;
Bc1 = sys_c.Bc1;
Bc2 = sys_c.Bc2;
Cc = sys_c.Cc;
Dc1 = sys_c.Dc1;
Dc2 = sys_c.Dc2;

% Intermediate matrices.
Dc1Dp = Dc1*Dp;
Z = eye(size(Dc1Dp,1)) - Dc1Dp;
DpiZDc1 = Dp*inv(Z)*Dc1;

if (all(abs(sys_c.Ac(:)) < 1e-6) && all(abs(sys_c.Bc1(:)) < 1e-6) && ...
        all(abs(sys_c.Bc2(:)) < 1e-6))
    % No controller state defined.
    
    % Closed loop A matrix.
    A = Ap + Bp*inv(Z)*Dc1*Cp;

    % Closed loop B matrix.
    B = Bp;

    % Closed loop C matrix.
    C = Cp;

    % Closed loop D matrix.
    D = Dp;
else
    % Controller state defined.

    % Closed loop A matrix.
    
    A = [Ap + Bp*inv(Z)*Dc1*Cp, Bp*inv(Z)*Cc;
         Bc1*(eye(size(DpiZDc1,1))+DpiZDc1)*Cp, Ac+Bc1*Dp*inv(Z)*Cc];
    % Closed loop B matrix.
    B = [Bp*inv(Z)*Dc2; Bc2 + Bc1*Dp*inv(Z)*Dc2];

    % Closed loop C matrix.
    C = [(eye(size(DpiZDc1,1))+DpiZDc1)*Cp, Dp*inv(Z)*Cc];

    % Closed loop D matrix.
    D = Dp*inv(Z)*Dc2;
end
