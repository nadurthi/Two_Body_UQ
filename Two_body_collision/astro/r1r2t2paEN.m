%  [p,a,dE] = r1r2t2paEN(r1v,r2v,tx,mu,shortorlong)
%
%   uses Newton's method
%
%			r1v, r2v are the two vectors
%			tx is the time of flight
%			mu is the two-body gravitational constant
%			shortorlong = 0|1 (short|long)
%		

function [p,a,dE] = r1r2t2paEN(r1v,r2v,tx,mu,shortorlong)
	r1   = norm(r1v);
	r2   = norm(r2v);
	dnu  = acos( (r1v'*r2v)/r1/r2 )
	if nargin == 5,
		if shortorlong == 1,
			dnu = 2*pi-dnu;
		end;
	end;
	k    = r1*r2*(1-cos(dnu))
	l    = r1+r2
	m    = r1*r2*(1+cos(dnu))
	ppi  = k/(l+sqrt(2*m))
	ppii = k/(l-sqrt(2*m))
	pn   = (ppi +   ppii)/2

	a   = m*k*pn/( (2*m-l*l)*pn^2 + 2*k*l*pn - k^2)
	f   = 1 - r2/pn *(1-cos(dnu))
	dE  = acos(1-r1/a*(1-f))
	g   = r1*r2*sin(dnu)/sqrt(mu*pn)
	t   = g + sqrt(a^3/mu)*(dE-sin(dE))

	dt   = t - tx
    its = 0;  MAXITS = 100;
	while ( abs(dt) > 1.0e-12 & its < MAXITS )
        its = its+1;
	    a   = m*k*pn/( (2*m-l*l)*pn^2 + 2*k*l*pn - k^2)
	    f   = 1 - r2/pn *(1-cos(dnu))
	    g   = r1*r2*sin(dnu)/sqrt(mu*pn)

        if ( a > 0 )
	        dE  = acos(1-r1/a*(1-f))
	        t   = g + sqrt(a^3/mu)*(dE-sin(dE))
            dtdp = -g/2/pn-3*a/2*(t-g)*(k^2+(2*m-l*l)*pn*pn)/m/k/pn^2+...
                   sqrt(a^3/mu)*2*k*sin(dE)/pn/(k-l*pn)
        else
            dF = acosh(1-r1*(1-f)/a)
            t = g + sqrt((-a)^3/mu)*(sinh(dF)-dF)
            dtdp = -g/2/pn-3*a/2*(t-g)*(k^2+(2*m-l*l)*pn*pn)/m/k/pn^2-...
                   sqrt((-a)^3/mu)*2*k*sinh(dF)/pn/(k-l*pn)
        end
        dp = (tx - t)/dtdp 
        pn = pn + dp
        dt = t-tx 
		%dtshow   = [dp t  t-tx ];
        pause
	end;
    show = [dp dt ];
    p = pn
    if (a > 0 )
	    dE = dE;
    else
        dE = dF;
    end
