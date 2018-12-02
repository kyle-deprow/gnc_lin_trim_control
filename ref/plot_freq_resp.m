function plot_freq_resp(freq_resp, w)
%**************************************************************************
%
% Generates various figures for the specified frequency response.
%
%
% Inputs:
%   freq_resp   Structure, frequency response data.  Refer to the file
%               "freq_analysis.m" for field definitions
%   w           n-element vector, frequency vector at which frequency
%               response data was collected (rps)
%
%
% Outputs:
%   MATLAB figures
%
%
% Version history:
%   2013.12.30  KRB  Initial release
%   2016.01.11  KRB  Placed Bode mag and phase on same figure
%
%**************************************************************************

% Constant
rtd = 180/pi;

% Extract data from frequency response structure.
L_in    = freq_resp.input.Lin;
wc      = freq_resp.input.LoopGainXover_rps;
RD_in   = freq_resp.input.ReturnDiff;
rd_min  = min(abs(RD_in));
SR_in   = freq_resp.input.StabRobust;
sr_min  = min(abs(SR_in));
sens    = freq_resp.output.S;
S_max   = max(abs(sens));
compsens= freq_resp.output.T;
T_max   = max(abs(compsens));
Gnois   = freq_resp.input.Gnoise;

% Draw circles.
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);
xx2=rd_min*cos(dd)-1;yy2=rd_min*sin(dd);

% Nyquist plot
figure,
plot(xx1,yy1,'k',real(L_in),imag(L_in),'b',xx2,yy2,'r-','LineWidth',2);grid
axis([-2 2 -2 2]);
text(.5,.5,['radius = ' num2str(rd_min)])
xlabel('Re(L)')
ylabel('Im(L)')
title('Nyquist')

xw_plot = 2*min(w);

% Loop Gain Bode magnitude and phase.
figure
subplot(2,1,1);
semilogx(w,20*log10(abs(L_in)),'LineWidth',2);grid
text(xw_plot,10.,['LGCF = ' num2str(wc)])
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title('Bode');
subplot(2,1,2);
semilogx(w,rtd*angle(L_in),'LineWidth',2);grid
xlabel('Frequency (rps)')
ylabel('Phase deg')

% Plot Return Difference 1+L magnitude.
figure
semilogx(w,20*log10(abs(RD_in)),'LineWidth',2);grid
text(xw_plot,10.,['min = ' num2str(rd_min)])
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title('|I+L| at input')

% Plot Stability Robustness 1+inv(L) magnitude.
figure
semilogx(w,20*log10(abs(SR_in)),'LineWidth',2);grid
text(xw_plot,10.,['min = ' num2str(sr_min)])
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title(' |I+inv(L)| at Plant Input')

% Plot sensitivity.
figure
semilogx(w,20*log10(sens),'LineWidth',2);grid
text(xw_plot,-10.,['max = ' num2str(S_max)])
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title(' |S| at Plant Output')

% Plot complementary sensitivity.
figure
semilogx(w,20*log10(compsens),'LineWidth',2);grid
text(xw_plot,-10.,['max = ' num2str(T_max)])
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title(' |T| at Plant Output')

% Plot noise to control magnitude.
figure
semilogx(w,20*log10(Gnois),'LineWidth',2);grid
xlabel('Frequency (rps)')
ylabel('Magnitude dB')
title('Noise-to-Control Transfer Function Matrix')