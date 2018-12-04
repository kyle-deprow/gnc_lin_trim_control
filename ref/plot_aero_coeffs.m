function plot_aero_coeffs(alphav_deg)

% Get the number of breakpoints.
n_aoa = numel(alphav_deg);

% Straight flight, various controls.
betav_deg = 0;
delev_deg = -20:5:20;
n_dele    = numel(delev_deg);

% Legend string.
leg_str = [];
for i = 1:n_dele
    str = sprintf('del = %d deg', delev_deg(i));
    leg_str = [leg_str, {str}];
end

% Loop over each Altitude.
for i_beta = 1:numel(betav_deg)
    % Get the beta.
    beta_deg = betav_deg(i_beta);
    
    % Instatiate a figure window
    figure
    set(gcf,'papertype','usletter')
    set(gcf,'paperposition',[0 0 8.5 11.0])
    set(gcf,'position',get(0,'ScreenSize'));
    orient landscape
    
    % Zero out memory for the data.
    zero_mat = zeros(n_aoa, n_dele);
    Cx = zero_mat;
    Cy = zero_mat;
    Cz = zero_mat;
    Cl = zero_mat;
    Cm = zero_mat;
    Cn = zero_mat;   
    
    % Loop over Speed.
    for i_def = 1:n_dele
        % Get the deflection.
        dele_deg = delev_deg(i_def);
        
        dela_deg = dele_deg;
        delr_deg = dele_deg;
        
        % Loop over Angle of Attack.
        for i_aoa = 1:n_aoa
            % Get the angle of attack.
            alpha_deg = alphav_deg(i_aoa);
            
            % Get Cx coefficient.
            Cx(i_aoa,i_def) = f16_cx(alpha_deg, ...
                                     dele_deg);       
            
            % Get Cy coefficient.
            Cy(i_aoa,i_def) = f16_cy(beta_deg, ...
                                     dela_deg, ...
                                     delr_deg);
            
            % Get Cz coefficient.
            Cz(i_aoa,i_def) = f16_cz(alpha_deg, ...
                                     beta_deg, ...
                                     dele_deg);
            
            % Get Cl coefficient.
            Cl(i_aoa,i_def) = f16_cl(alpha_deg, beta_deg) + ...
                              f16_cl_dail(alpha_deg, beta_deg)*(dela_deg/20) + ...
                              f16_cl_drud(alpha_deg, beta_deg)*(delr_deg/30);
            
            % Get Cm coefficient.
            Cm(i_aoa,i_def) = f16_cm(alpha_deg, ...
                                     dele_deg);
            
            % Get Cn coefficient.
            Cn(i_aoa,i_def) = f16_cn(alpha_deg, beta_deg) + ...
                              f16_cn_dail(alpha_deg, beta_deg)*(dela_deg/20) + ...
                              f16_cn_drud(alpha_deg, beta_deg)*(delr_deg/30);
        end
    end
    
    % Plot the data.
    subplot(3,2,1)
    plot(alphav_deg,Cx);grid on;
    xlabel('Angle of Attack, deg');
    ylabel('Cx');
    legend(leg_str,'Location','NorthEastOutside');

    
    subplot(3,2,2)
    plot(alphav_deg,Cl);grid on;
    xlabel('Angle of Attack, deg');
    ylabel('Cl');
    legend(leg_str,'Location','NorthEastOutside');
    
    subplot(3,2,3)
    plot(alphav_deg,Cy);grid on;
    xlabel('Angle of Attack, deg');
    ylabel('Cy');
    legend(leg_str,'Location','NorthEastOutside');
    
    subplot(3,2,4)
    plot(alphav_deg,Cm);grid on;
    xlabel('Angle of Attack, deg');
    ylabel('Cm');
    legend(leg_str,'Location','NorthEastOutside');
    
    subplot(3,2,5)
    plot(alphav_deg,Cz);grid on;
    xlabel('Angle of Attack, deg');
    ylabel('Cz');
    legend(leg_str,'Location','NorthEastOutside');
    
    subplot(3,2,6)
    plot(alphav_deg,Cn);grid on;
    xlabel('Angle of Attack, deg');
    ylabel('Cn');
    legend(leg_str,'Location','NorthEastOutside');
end
