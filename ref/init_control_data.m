function control_data = init_control_data()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function initializes a data structure for n number of design points.
%
%
% Inputs:
%   n       scalar, number of data structures to initialize
%
%
% Outputs:
%   data_point  structure array, array of data
%
%
% Version history:
%   2016.10.08  KRB  Initial release
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define fields for transient response.
control_data.PO     = [];
control_data.tr_sec = [];
control_data.ts_sec = [];

% Define fields for control.
control_data.Q = [];
control_data.R = [];
control_data.maxu = [];

% Define fields for frequency response.
control_data.sr = [];
control_data.rd = [];
control_data.wc_Hz = [];
