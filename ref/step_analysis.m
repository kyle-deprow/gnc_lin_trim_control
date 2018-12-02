function s = step_analysis(sys, i_in, i_out, step_mag)
%**************************************************************************
%
% This function computes step response statistics for a system.  This
% relates a user specified input to a set of user specified outputs.
%
%
% Inputs:
%   sys          structure, state space structure of the system
%                   .a,     A matrix
%                   .b,     B matrix
%                   .c,     C matrix
%                   .d,     D matrix
%   i_in        n_in element vector, indices of the system input space to
%               compute transient responses against
%   i_out       n_out element vector, indices of the system output space 
%               for which compute transient responses
%   step_mag    scalar, step input magnitude
%
%
% Outputs:
%   s       n_in*n_out instance structure, array containing the transient
%           response data for each input/output pairs
%               .input,     input channel
%               .output,    output channel
%               .tr63_sec,  63% rise time (sec)
%               .tr90_sec,  90% rise time (sec)
%               .tp_sec,    time of first peak (sec)
%               .PO,        percent overshoot
%               .tu_sec,    time of first undershoot (sec)
%               .PU,        percent undershoot
%               .ts98_sec,  98% settling time (sec)
%
%
% Version history:
%   2015.12.22  KRB  Initial release
%
%**************************************************************************

% Get the number of inputs and outputs to process.
n_in = numel(i_in);
n_out= numel(i_out);

% Extract the system matrices.
a = sys.a;
b = sys.b;
c = sys.c;
d = sys.d;

% Initialize the structure.
for i = 1:n_in*n_out
    s(i).input = [];
    s(i).output = [];
    s(i).tr63_sec = [];
    s(i).tr90_sec = [];
    s(i).tp_sec = [];
    s(i).PO = [];
    s(i).tu_sec = [];
    s(i).PU = [];
    s(i).ts98_sec = [];
end    

% Loop through the input / output combinations.
icnt = 0;
for i = 1:n_in
    % Get the input index.
    i_index = i_in(i);
    
    for j = 1:n_out
        % Get the output index.
        o_index = i_out(j);
        
        % Build the matrices for analysis.
        B = b(:,i_index);
        C = c(o_index, :);
        D = d(o_index, i_index);
        
        % Build the system.
        sys_step = ss(a,B,C,D);
        
        % Perform the step response.
        [y, time_s] = step(sys_step*step_mag);
        
        % Build a time vector for analysis.
        t_s = 0:0.001:time_s(end);
        
        % Perform the step response.
        y = step(sys_step*step_mag, t_s);
        
        % Get the steady state value.  Assuming the steady state value is
        % at the end of the simulation.
        yss = y(end);
        
        % Compute the 63% rise time (or set to 0).
        i_y63 = find(y >= 0.63*yss,1);
        if isempty(i_y63)
            % Response did not get to the 63% level
            tr63 = 0;
        else
            % Get the time associated with the 63% level.
            tr63 = t_s(i_y63);
        end
        
        % Compute the 90% rise time (or set to 0).
        i_y10 = find(y >= 0.10*yss,1);
        i_y90 = find(y >= 0.90*yss,1);
        if isempty(i_y90)
            % Response did not get to the 90% level.
            tr90 = 0;
        else
            % Compute the 90% rise time between the 10% to 90% steady state
            % value.
            tr90 = t_s(i_y90) - t_s(i_y10);
        end
        
        % Compute Percent Overshoot (or set to 0).
        ymax = max(y);
        if ymax > yss
            % Overshoot occurred.  Compute the peak time.
            tp = t_s(find(y == ymax, 1));
            
            % Compute the Percent overshoot.
            PO = 100*((ymax-yss)/yss);
        else
            % Overshoot did not occur.
            tp = 0;
            PO = 0;
        end
        
        % Compute Percent Undershoot (or set to 0).
        ymin = min(y);
        if ymin < 0
            % Undershoot occurred.  Compute the undershoot time.
            tu = t_s(find(y == ymin, 1));
            
            % Get the y value at undershoot.
            US = 100*(abs(ymin)/yss);
        else
            % Undershoot did not occur.
            tu = 0;
            US = 0;
        end
        
        % Compute the 98% settling time (or set to 0).
        ss_error = abs(y - yss);
        i_e = find(ss_error > 0.02*yss, 1, 'last');
        if i_e < numel(y)
            % Settling occurred.  Compute the settling time.
            ts98 = t_s(i_e) + ...
                   (t_s(i_e+1)-t_s(i_e))*(0.02*yss-ss_error(i_e)) / ...
                   (ss_error(i_e+1)-ss_error(i_e));
        else
            % Settling did not occur.
            ts98 = 0;
        end
        
        % Store the data.
        icnt = icnt + 1;
        s(icnt).input = i_index;
        s(icnt).output = o_index;
        s(icnt).tr63_sec = tr63;
        s(icnt).tr90_sec = tr90;
        s(icnt).tp_sec = tp;
        s(icnt).PO = PO;
        s(icnt).tu_sec = tu;
        s(icnt).PU = US;
        s(icnt).ts98_sec = ts98;
    end
end