%  t= = tofr1r2p(r1v,r2v,p,mu,shortorlong)
%
%   one iteration of p-iteration method
%
%			r1v, r2v are the two vectors
%			mu is the two-body gravitational constant
%			shortorlong = 0|1 (short|long)
%		

function t = tofr1r2p(r1v,r2v,p,mu,shortorlong)
	r1   = norm(r1v)
	r2   = norm(r2v)
   dnu  = acos( (r1v'*r2v)/r1/r2 )
   if nargin == 5,
		if shortorlong == 1,
			dnu = 2*pi-dnu
		end;
	end;
	k    = r1*r2*(1-cos(dnu))
	l    = r1+r2
	m    = r1*r2*(1+cos(dnu))
   
	a    = m*k*p /( (2*m-l*l)*p^2 + 2*k*l*p - k^2)
   f    = 1 - r2/p *(1-cos(dnu))
   fdot = sqrt(mu/p)*tan(dnu/2)*((1-cos(dnu))/p-1/r1-1/r2)
   dE   = atan2(-r1*r2*fdot/sqrt(mu*a),1-r1/a*(1-f))
   if (dE < 0), dE = dE+2*pi, end
	g    = r1*r2*sin(dnu)/sqrt(mu*p)
   t    = g + sqrt(a^3/mu)*(dE-sin(dE))
   dE-sin(dE)
   TP = 2*pi*sqrt(a^3/mu)
