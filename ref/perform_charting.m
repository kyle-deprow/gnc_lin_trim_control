function cont_data = perform_charting(sys_p, Aw, Bw, Creg, Qw, Rw, q, w)
% AENG - 555 : Guidance and Control of Aerospace Vehicles
% Saint Louis University
% Ken Buckholtz, Ph.D.
%
% Version history:
%   November 2, 2017   Initial release
%   November 21, 2018  Add "try-catch" in case system cannot be stabilized

% Initialize the output data structure.
cont_data = init_control_data();

% Loop through the scale factors and collect up performance metrics.
PO = zeros(numel(q),1);
tr = PO;
ts = PO;
wc = PO;
maxu=PO;
srmin = PO;
rdmin = PO;

% Number of channels to regulate.
nreg = size(Creg,1);

% Get the dimensions of the problem.
[n,m] = size(sys_p.Bp);

% Define magnitude of step.
step_mag = 0.1;

for i = 1:numel(q)
    % Set the integral error state weights to the current instance of the
    % scale factor.
    for j = 1:nreg
        Qw(j,j) = q(i);
    end
    
    % Solve the LQR problem.
    try 
        [Kw, S, E] = lqr(Aw, Bw, Qw, Rw);
    catch
        PO(i) = NaN;
        tr(i) = NaN;
        ts(i) = NaN;
        wc(i) = NaN;
        rdmin(i) = NaN;
        srmin(i) = NaN;
        maxu(i) = NaN;
        continue
    end
    
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
    
    % Get the eigenvalues.
    ee = eig(Acl);
    if any(real(ee) > 0)
        continue
    end
    
    % Perform frequency response from control (channel 1) to the ball
    % position output (output 1).
    freq_response = freq_analysis(sys_p, sys_c, w, 1, 3);
    
    % Get the Loop Gain crossover frequency.
    wc(i) = freq_response.input.LoopGainXover_rps;
    
    % Get the Return Difference, I + L, and its min magnitude.
    rd = freq_response.input.ReturnDiff;
    rdmin(i) = min(abs(rd));
    
    % Get the Stability Robustness, I + inv(L), and its min magnitdue.
    sr = freq_response.input.StabRobust;
    srmin(i) = min(abs(sr));
    
    % Perform a step response.
    syscl = ss(Acl,Bcl,Ccl,Dcl);
    [y, time_sec, x] = step(syscl*step_mag);
    
    % Reconstruct the control and pick the biggest magnitude.
    u = -Kw(:,[nreg+1:nreg+n,1:nreg])*x';
    maxu(i) = max(abs(u));
    
    if cond(Acl) > 50000
        % Compute step information for the system.  Don't use stepinfo() as
        % it does not deliver correct results for multi-outputs.
        s = step_analysis(syscl, 1, 3, step_mag);
        
        % Get the transient metrics.
        PO(i) = s.PO;
        tr(i) = s.tr90_sec;
        ts(i) = s.ts98_sec;
    else
        sout = stepinfo(syscl*step_mag);
        
        % Get the transient metrics.
        PO(i) = sout(3).Overshoot;
        tr(i) = sout(3).RiseTime;
        ts(i) = sout(3).SettlingTime;
    end
end

% Save control data.
cont_data.tr_sec = tr;
cont_data.ts_sec = ts;
cont_data.PO = PO;
cont_data.sr = srmin;
cont_data.rd = rdmin;
cont_data.maxu = maxu;
cont_data.wc_Hz = wc/(2*pi);


% Generate LQR design charts.
figure
subplot(3,2,1);
semilogx(q, PO); grid on;
xlabel('q');
ylabel('PO');
subplot(3,2,2);
semilogx(q, tr); grid on;
xlabel('q');
ylabel('t_r, sec');
subplot(3,2,3);
semilogx(q, ts); grid on;
xlabel('q');
ylabel('t_s, sec');
subplot(3,2,4);
semilogx(q, wc/(2*pi)); grid on;
xlabel('q');
ylabel('f_c, Hz');
subplot(3,2,5);
semilogx(q,maxu); grid on;
xlabel('q');
ylabel('max u');
subplot(3,2,6);
semilogx(q,rdmin,'b-',q,srmin,'b--'); grid on;
xlabel('q');
ylabel('minRD / minSR');
