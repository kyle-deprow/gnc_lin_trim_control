function plot_trim_results(design_point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates plots of the trim results
%
%
% Inputs:
%   design_point    n-element structure array, flight condition data
%
%
% Outputs:
%   plots
%
%
% Version history:
%   2018.11.04  KRB  Initial release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of data points.
n_pts = numel(design_point);

% Create a vector of flight conditions.
idx_fc_sel = 1:n_pts;

% Create a vector of flight conditions that trimmed.
idx_valid_trim = find([design_point.trim_flag] == 1);

% define indices of invalid trim solutions
idx_invalid_trim = find([design_point.trim_flag] ~= 1);

% define subset of selected flight conditions with valid trim conditions
idx_to_plot = intersect(idx_valid_trim,idx_fc_sel);

% define subset of selected flight conditions with invalid trim conditions.
idx_invalid_to_plot = intersect(idx_invalid_trim, idx_fc_sel);

% define independent variable for plots.  Creates a set of indices for
% those points that trimmed and a set for those that did not trim.
[~,~,indep_var] = intersect(idx_to_plot,idx_fc_sel);
[~,~,indep_var_invalid] = intersect(idx_invalid_to_plot,idx_fc_sel);
xlabel_str = 'flight cond #';


% define plot options
linestyle = 'bx';
invalid_linestyle = 'ro';
linewidth  = 2;

% set x-axis limits
xlim_min = min(idx_fc_sel);
xlim_max = max(idx_fc_sel);
if xlim_min == xlim_max
    xlim_min = xlim_min - 1;
    xlim_max = xlim_max + 1;
end

% Get the trim results.
xtrim = [design_point.xtrim];
utrim = [design_point.utrim];


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 1: alt/mach/alpha/fuel/Nz
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure
set(gcf,'papertype','usletter')
set(gcf,'paperposition',[0 0 8.5 11.0])
orient landscape

% Plot Altitude
subplot(4,1,1)
hold on
scatter(indep_var,[design_point(idx_to_plot).alt_ft],...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,[design_point(idx_invalid_to_plot).alt_ft],...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('alt (ft)'), xlim([xlim_min xlim_max]), grid on
title('Trim Results: Trimmed (blue), Did not Trim (red)');

% Plot Speed
subplot(4,1,2)
hold on
scatter(indep_var,[design_point(idx_to_plot).V_fps],...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,[design_point(idx_invalid_to_plot).V_fps],...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('Airspeed (ft/sec)'), xlim([xlim_min xlim_max]), grid on

% Plot Angle of Attack
subplot(4,1,3)
hold on
scatter(indep_var,[design_point(idx_to_plot).alpha_rad]*(180/pi),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,[design_point(idx_invalid_to_plot).alpha_rad]*(180/pi),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('Angle of Attack (deg)'), xlim([xlim_min xlim_max]), grid on

% Plot Normal accel
subplot(4,1,4)
hold on
scatter(indep_var,-[design_point(idx_to_plot).Az_g],...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,-[design_point(idx_invalid_to_plot).Az_g],...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('Nz (+g up)'), xlim([xlim_min xlim_max]), grid on
xlabel(xlabel_str), linkaxes(datachildren(gcf),'x')


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 2: longitudinal control
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure
set(gcf,'papertype','usletter')
set(gcf,'paperposition',[0 0 8.5 11.0])
orient landscape

% Plot 1-g normal acceleration line and normal acceleration
subplot(4,1,1)
plot([indep_var(1) indep_var(end)],[1,1],'r');
hold on
scatter(indep_var,-[design_point(idx_to_plot).Az_g],...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,-[design_point(idx_invalid_to_plot).Az_g],...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('NzCG (g+up)'), xlim([xlim_min xlim_max]), grid on
title('Trim Results: Trimmed (blue), Did not Trim (red)');

% Plot Angle of Attack
subplot(4,1,2)
hold on
scatter(indep_var,[design_point(idx_to_plot).alpha_rad]*(180/pi),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,[design_point(idx_invalid_to_plot).alpha_rad]*(180/pi),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('alpha (deg)'), xlim([xlim_min xlim_max]), grid on

% Plot Elevator deflection
subplot(4,1,3)
hold on
dele_deg = utrim(2,:);
scatter(indep_var,dele_deg(idx_to_plot),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,dele_deg(idx_invalid_to_plot),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('dele (deg)'), xlim([xlim_min xlim_max]), grid on

% Plot trim thrust
subplot(4,1,4)
hold on
delt = utrim(1,:);
scatter(indep_var,delt(idx_to_plot),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,delt(idx_invalid_to_plot),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('delt'), xlim([xlim_min xlim_max]), grid on
xlabel(xlabel_str), linkaxes(datachildren(gcf),'x')


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 3: lateral-directional control
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure
set(gcf,'papertype','usletter')
set(gcf,'paperposition',[0 0 8.5 11.0])
orient landscape

% Plot Angle of Sideslip
subplot(3,1,1)
hold on
beta_deg = xtrim(3,:)*(180/pi);
scatter(indep_var,beta_deg(idx_to_plot),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,beta_deg(idx_invalid_to_plot),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('beta (deg)'),
xlim([xlim_min xlim_max]), grid on
title('Trim Results: Trimmed (blue), Did not Trim (red)');

% Plot Aileron deflection
subplot(3,1,2)
hold on
dela_deg = utrim(3,:);
scatter(indep_var,dela_deg(idx_to_plot),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,dela_deg(idx_invalid_to_plot),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('dela (deg)'), xlim([xlim_min xlim_max]), grid on

% Plot Rudder deflection
subplot(3,1,3)
hold on
delr_deg = utrim(4,:);
scatter(indep_var,delr_deg(idx_to_plot),...
        linestyle,'LineWidth',linewidth)
scatter(indep_var_invalid,delr_deg(idx_invalid_to_plot),...
        invalid_linestyle,'LineWidth',linewidth)
hold off
ylabel('delr (deg)'), xlim([xlim_min xlim_max]), grid on
xlabel(xlabel_str), linkaxes(datachildren(gcf),'x')

