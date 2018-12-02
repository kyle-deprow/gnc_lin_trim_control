function [a, b] = linearize_nonlinear_model(name, x, u, u1)

% Get the size of the state and controls.
n = numel(x);
m = numel(u);
tol=1e-6; 
time=0.;

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
   for i=1:60
      xt(j)=x(j)+dx(j);
      xd1= feval(name,time,xt,[u;u1]);
      xt(j)=x(j)-dx(j);
      xd2= feval(name,time,xt,[u;u1]);
      a(:,j)= (xd1-xd2)'/(2*dx(j));
      if max( abs(a(:,j)-last)./abs( a(:,j) + 1e-12 ) )<tol;
         break
      end
      dx(j)= 0.125*dx(j);
      last = a(:,j);
   end
   %column=j
   iteration=i;
   if iteration==20
      disp(['not converged on A, column',num2str(j)])
   end
end


% Create the Jacobian for the B matrix by numerically computing the partial
% derivatives of the state equations with respect to each input.
du=0.1*u;
for i=1:m
   if du(i)==0.0;
      du(i)=0.1;
   end
end
last=zeros(n,1); b=zeros(n,m);
for j=1:m
   last=zeros(n,1);
   usave=u;
   for i=1:10
      u(j)=usave(j)+du(j);
      xd1= feval(name,time,x,[u;u1]);
      u(j)=usave(j)-du(j);
      xd2= feval(name,time,x,[u;u1]);
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
