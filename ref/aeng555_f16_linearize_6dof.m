% Saint Louis University
% Ken Buckholtz, Ph.D.
%
% Version history:
%   September 21, 2011  Initial release
%   November  13, 2014  Add the computation of the B matrix for Lat/Dir
%
% This script linearizes the nonlinear F-16 model at the prescribed flight
% condition.
%
% Data taken from :
%   Aircraft Control and Simulation, 2nd edition
%   by Stevens and Lewis
set(0,'DefaultTextInterpreter','tex');
clear
close all


% Nonlinear dynamics model name.
tol=1e-6; 
time=0.;
name = 'f16_nonlinear_model';
% Set the Xcg percentage of cbar.
xcg = 0.3;
% Trim condition on the state vector.  Obtain this from the trim routine.
x = [500
     0.0349
         0
         0
         0
         0
         0
         0
         0
         0
         0
         1
    6.4803];
xretain = [1:9,13];

% Trim condition on the control vector.  Obtain this from the trim routine.
u = [0.0997888
    -1.866
    0.0
    0.0];

% Set the number of states and controls.
n = numel(x);
mm= numel(u);
m = mm;

% Set the number of outputs.
p = 3;

%{
% Create the Jacobian for the A matrix by numerically computing the partial
% derivatives of the state equations with respect to each state.
dx=0.1*x;
for i=1:n
   if dx(i)==0.0;
      dx(i)=0.1;
   end
end
last=zeros(n,1);  a=zeros(n,n);
for j=1:n
   xt=x;
   for i=1:10
      xt(j)=x(j)+dx(j);
      xd1= feval(name,time,xt,[u;xcg]);
      xt(j)=x(j)-dx(j);
      xd2= feval(name,time,xt,[u;xcg]);
      a(:,j)= (xd1-xd2)'/(2*dx(j));
      if max( abs(a(:,j)-last)./abs( a(:,j) + 1e-12 ) )<tol;
         break
      end
      dx(j)= 0.5*dx(j);
      last = a(:,j);
   end
   %column=j
   iteration=i;
   if iteration==10
      disp(['not converged on A, column',num2str(j)])
   end
end


% Create the Jacobian for the B matrix by numerically computing the partial
% derivatives of the state equations with respect to each input.
du=0.1*u;
for i=1:mm
   if du(i)==0.0;
      du(i)=0.1;
   end
end
last=zeros(n,1); b=zeros(n,m);
for j=1:mm
   usave=u;
   for i=1:10
      u(j)=usave(j)+du(j);
      xd1= feval(name,time,x,[u;xcg]);
      u(j)=usave(j)-du(j);
      xd2= feval(name,time,x,[u;xcg]);
      b(:,j)= (xd1-xd2)'/(2*du(j));
      if max( abs(b(:,j)-last)./abs( b(:,j) + 1e-12 ) )<tol;
         break
      end
      du(j)= 0.5*du(j);
      last = b(:,j);
   end
   %column=j
   iteration=i;
   if iteration==10
      disp('not converged on B, column',j)
   end
end
%}
% Obtain the linearized state equation (A,B).
[a, b] = linearize_nonlinear_model(name, x, u, xcg);

% Create A and B matrix for those states retained.
A = a(xretain,xretain);
B = b(xretain,:);
% Create the Jacobian for the C matrix by numerically conputing the partial
% derivative of the output with respect to each state.
dx=0.1*x;
for i=1:n  
   if dx(i)==0.0;
      dx(i)=0.1;
   end
end
last=zeros(p,1);  c=zeros(p,n);
%xd1=zeros(n);  xd2=zeros(n);
for j=1:n
   xt=x;
   for i=1:10
      xt(j)=x(j)+dx(j);
      [xd1, az1] = feval(name,time,xt,[u;xcg]);
      xt(j)=x(j)-dx(j);
      [xd2, az2] = feval(name,time,xt,[u;xcg]);
      c(1,j) = (az1-az2)/(2*dx(j));
      if max( abs(c(:,j)-last)./abs( c(:,j) + 1e-12 ) )<tol;
         break
      end
      dx(j)= 0.5*dx(j);
      last = c(:,j);
   end
   %column=j
   iteration=i;
   if iteration==10
      disp(['not converged on C, column',num2str(j)])
   end
end
% Insert output for q state (deg/sec), and AOA state (deg).

c(2,8) = 180/pi;
c(3,2) = 180/pi;
% Compute the Jacobian of the D matrix by numerically computing the partial
% derivative of the output with respect to each input.
du=0.1*u;
for i=1:mm
   if du(i)==0.0;
      du(i)=0.1;
   end
end
last=zeros(p,1); d=zeros(p,m);
for j=1:mm
   usave=u;
   for i=1:10
      u(j)=usave(j)+du(j);
      [xd1,az1,ay1] = feval(name,time,x,[u;xcg]);
      u(j)=usave(j)-du(j);
      [xd2,az2,ay2] = feval(name,time,x,[u;xcg]);
      d(1,j) = (az1-az2)/(2*du(j));
      if max( abs(d(:,j)-last)./abs( d(:,j) + 1e-12 ) )<tol;
         break
      end
      du(j)= 0.5*du(j);
      last = d(:,j);
   end
   %column=j
   iteration=i;
   if iteration==10
      disp('not converged on D, column',j)
   end
end
% Create C and D matrix for those states retained.
C = c(:,xretain);
D = d;
% Create the Longitudinal Modes Matrix, A, by only using the states
% associated with Vt, AOA, theta and q.
Am = A([1 2 5 8],[1 2 5 8]);

% Create the Longitudinal Input Matrix, B, by only using the states
% associated with Vt, AOA, theta and q, and controls delT, dele.
Bm = B([1 2 5 8], 1:2);

% Compute the eigenvalues and eigenvectors.
[V,Eig] = eig(Am);

% Analysis of transfer function between pitch rate, q, output and elevator,
% de, input.
B([1,2,5,8],2)
C(2,[1,2,5,8])
D(2,2)
[num,den] = ss2tf(Am,B([1 2 5 8],2),C(2,[1 2 5 8]),D(2,2));
G = tf(num,den);
G = zpk(G);
% Short period approximation.
numa = -11.463*[1 1.034];
dena = [1 2.516 4.023];
Ga = tf(numa,dena);
Ga = zpk(Ga);
% Bode plots.
bode(G,'b-',Ga,'b--');grid on
title('Bode Plots for q to \delta_e for Short + Phugoid (solid) and Short Approx(dashed)');


% Create the Lateral/Directional Modes Matrix, A, by only using the states
% associated with AOS, phi, p and r.
Ald = A([3 4 7 9], [3 4 7 9]);

% Create the Lateral/Directional Input Matrix, B, by only using the states
% associated with AOS, phi, p and r, and controls, dela, delr.
Bld = B([3 4 7 9],3:4);

% Compute the eigenvalues and eigenvectors.
[Vld, Eigld] = eig(Ald);

