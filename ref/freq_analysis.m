function [freq_data] = freq_analysis(sys_p, sys_c, w, i_in, i_out)
%**************************************************************************
%
% Function performs a frequency analysis for a MIMO system for a particular
% input / output pair.  The analysis is performed at both the plant input
% and the plant output.
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
%   freq_data   structure, frequency response data
%               Fields:
%               .input      structure, frequency response broken at plant
%                           input
%               .output     structure, frequency response broken at plant
%                           output
%               .system     structure, data associated with closed loop
%                           system of plant and controller
%
%
% Version history:
%   2013.12.30  KRB  Initial release
%   2014.03.17  KRB  Fixed Cin in computing Lin 
%
%**************************************************************************

% Get the number of frequencies.
nf = numel(w);

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

% Initialize input frequency analysis substructure.
freq_data.input.Lin = zeros(nf,1);
freq_data.input.ReturnDiff= zeros(nf,1);
freq_data.input.StabRobust= zeros(nf,1);
freq_data.input.LoopGainXover_rps= 0;
freq_data.input.GMRD = [0,0];
freq_data.input.PMRD = [0,0];
freq_data.input.GMSR = [0,0];
freq_data.input.PMSR = [0,0];
freq_data.input.Gnoise = zeros(nf,1);

% Initialize output frequency analysis substructure.
freq_data.output.S = zeros(nf,1);
freq_data.output.T = zeros(nf,1);

% Initialize system substructure.
freq_data.system.OLEvalues = [];
freq_data.system.CLEvalues = [];

% Close the loop
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [(Ap+Bp*Z*Dc1*Cp)               (Bp*Z*Cc);
       (Bc1*Cp+Bc1*Dp*Z*Dc1*Cp) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [(Bp*Z*Dc2);
     (Bc2+Bc1*Dp*Z*Dc2)];
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);


%%%%% EIGENVALUE ANALYSIS %%%%%

% Compute the eignevalues of the open loop plant, and closed loop system.
freq_data.system.OLEvalues = eig(Ap);
freq_data.system.CLEvalues = eig(Acl);


%%%%% LOOP GAIN AT PLANT INPUT ANALYSIS %%%%%

%Create SS model of loop gain at the plant input.
Ain = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac];
Bin = [ Bp; Bc1*Dp];
Cin = -[ Dc1*Cp Cc];
Din = -[ Dc1*Dp];

RD_in= zeros(nf,1);
SR_in= zeros(nf,1);
Gnois= zeros(nf,1);

for i = 1:nf
    % s = j*omega
    s = sqrt(-1)*w(i);

    % Controller object.
    KK = Cc*inv(s*eye(size(Ac))-Ac)*Bc1+Dc1;

    % Loop Gain.
    L_in(i) = Cin(i_in,:)*inv(s*eye(size(Ain))-Ain)*Bin(:,i_in)+Din(i_in,i_in);

    % Return difference.
    RD_in(i) = 1.+L_in(i);

    % Difference at actuator input.
    SR_in(i) = 1.+1./L_in(i);

    % Noise transfer function.
    Gnois(i)   = max(svd(KK));
end

% Compute crossover frequency.
magdb = 20.*log10(abs(L_in));
wc = FindCrossover(magdb,w);

% Compute the return difference and stability robustness MIMO margins.
rtd = 180/pi;
rdm = min(abs(RD_in));
srm = min(abs(SR_in));
rd_gm = [ (1/(1+rdm)) (1/(1-rdm)) ];
sr_gm = [ (1-rdm) (1+rdm) ];
rd_gm = 20*log10(rd_gm);
sr_gm = 20*log10(sr_gm);

rd_pm = rtd*2*asin(rdm/2);
sr_pm = rtd*2*asin(srm/2);

% Store.
freq_data.input.Lin = L_in;
freq_data.input.ReturnDiff= RD_in;
freq_data.input.StabRobust= SR_in;
freq_data.input.LoopGainXover_rps= wc;
freq_data.input.GMRD = rd_gm;
freq_data.input.PMRD = rd_pm;
freq_data.input.GMSR = sr_gm;
freq_data.input.PMSR = sr_pm;
freq_data.input.Gnoise = Gnois;


%%%%% LOOP GAIN AT PLANT OUTPUT ANALYSIS %%%%%

%SS model of loop gain at the plant output
Aout = [ Ap Bp*Cc;  0.*Bc1*Cp Ac];
Bout = [ Bp*Dc1; Bc1];
Cout = [ Cp Dp*Cc];
Dout = [ Dp*Dc1];

sens = zeros(nf,1);
compsens = zeros(nf,1);
for i = 1:nf
    % s = j*omega
    s = sqrt(-1)*w(i);
    
    % Loop at Plant output (Az).
    Lout = Cout(i_out,:)*inv(s*eye(size(Aout))-Aout)*Bout(:,i_in)+Dout(i_out,i_in);

    % Sensitivity
    Sout = inv(eye(size(Lout))+Lout);

    % Complementary sensitivity
    Tout = Sout*Lout;

    % Get max SV for each sensitivity.
    sens(i) = max(svd(Sout));
    compsens(i) = max(svd(Tout));
end

% Store.
freq_data.output.S = sens;
freq_data.output.T = compsens;




function t1 = FindCrossover(a,t)
%find the value of t where a crosses zero
n=numel(a);
j=0;
while j <= n,
  if a(n-j) > 0.,
    i=n-j;
    j=n+1;
  end;
  j=j+1;
end


pp=inv([t(i) 1.;t(i+1) 1.])*[a(i);a(i+1)];
t1=-pp(2)/pp(1);
