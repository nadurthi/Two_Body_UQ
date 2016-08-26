function [tVector,x] = rk4(deriv,tVector,x0)

% deriv = name of derivative m-file
% tVector = vector of time
% x0 = state, or initial conditions

% Note that this method requires a constant step-size
dt = tVector(2) - tVector(1);

x0 = x0(:);
x(1,:) = x0';

dt2 = dt/2;
dt6 = dt/6;

for i=2:length(tVector)
   
   t0 = tVector(i-1);
   t1 = tVector(i);
   
   for t = t0:dt:t1-dt
      tt = t + dt2;
      dxdt = feval(deriv,t,x0);
      x0t = x0 + dt2*dxdt;
      dx0t = feval(deriv,tt,x0t);
      x0t = x0 + dt2*dx0t;
      dxmold = feval(deriv,tt,x0t);
      x0t = x0 + dt*dxmold;
      dxm = dx0t + dxmold;
      dx0t = feval(deriv,t + dt,x0t);
      x0 = x0 + dt6*(dxdt + dx0t + 2*dxm);
   end
   
   x(i,:) = x0';
   
end
