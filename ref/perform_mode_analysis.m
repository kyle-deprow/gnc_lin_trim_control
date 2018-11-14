function [Msens, Mmetric] = perform_mode_analysis(A)
% Function computes the mode sensitivities for each state associated with
% the matrix A.  This gives an indication of which states influence each
% mode.  This function also computes a set of metrics for each mode.  Note
% that frequency compoents of a mode exist if eigenvectors have complex
% conjugate pairs.
%
%
% Inputs:
%   A       [nxn] matrix, System state matrix
%
%
% Outputs:
%   Msens   [nxn] matrix, Mode sensitivity matrix
%   Mmetric [3xn] matrix, Mode metrics.
%
%
% Version history:
%   2015.01.20  KRB  Initial release
%   2015.05.05  KRB  Added processing for purely oscillatory modes
%

% Get the dimension of the square system matrix.
n = size(A,1);

% Initialize the output.
Msens = zeros(n,n);
Mmetric = zeros(3,n)-999;

% Compute the eigenvectors, V = [v1 v2 ... vn], for the system.  Also
% returns an eigenvalue matrix, D, where the associated eigenvalues are
% down the main diagonal.
[V, D] = eig(A);

% Compute the inverse of matrix V.
Vinv = inv(V);

% Compute the sensitivities of each mode.
for i = 1:n
    % Build a diagonal matrix of the ith column of the inverse eigenvector
    % matrix.
    C = diag(Vinv(:,i));
    
    % Compute the modes of the ith state by multiplying the ith row of V by
    % the diagonal matrix of Vinv.  Take the magnitude in case the mode has
    % a frequency component (i.e. imaginary value).
    Msens(i,:) = abs(V(i,:)*C);
end


% Compute the metics of each mode.  A complex conjugate pair will realize
% the same metrics.
for i = 1:n
    % Get the real and imaginary components of the mode.
    Relambda = real(D(i,i));
    Imlambda = imag(D(i,i));
    
    % Check if mode is purely oscillatory.
    is_oscillation = abs(Relambda) < 1e-6;
    
    % Check if mode is purely real.
    is_real = abs(Imlambda) < 1e-6;
    
    % Compute the time to half amplitude for non-oscillatory modes.
    thalf = 0;
    if ~is_oscillation
        thalf = log(1/2) / Relambda;

        % Store the time to half amplitude.
        Mmetric(1,i) = thalf;
    end
    
    % If this eigenvalue is real, then there is no frequency component.
    if is_real
        continue
    end
    
    % Compute the period of the mode.
    T = (2*pi)/abs(Imlambda);
    
    % Store the period.
    Mmetric(2,i) = T;
    
    % Compute the ratio of time to half amplitude to the period.
    if ~is_oscillation
        Nhalf = thalf/T;

        % Store the ratio.
        Mmetric(3,i) = Nhalf;
    end
    Mmetric
end
