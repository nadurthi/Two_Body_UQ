%  [p,a,dE] = r1r2st2paE(r1,r2,dnu,tx,mu)
%   scalar version
%   uses the p-iteration method
%
%			r1, r2 are the two radii
%        dnu is the change in true anomaly
%			tx is the time of flight
%			mu is the two-body gravitational constant
%		

function [p,a,dE] = r1r2t2paE(r1,r2,dnu,tx,mu)
	k    = r1*r2*(1-cos(dnu));
	l    = r1+r2;
	m    = r1*r2*(1+cos(dnu));
	ppi  = k/(l+sqrt(2*m));
	ppii = k/(l-sqrt(2*m));
	po   = (2*ppi +   ppii)/3;
	p1   = (  ppi + 2*ppii)/3;

	ao   = m*k*po/( (2*m-l*l)*po^2 + 2*k*l*po - k^2);
	fo   = 1 - r2/po *(1-cos(dnu));
	dEo  = acos(1-r1/ao*(1-fo));
	go   = r1*r2*sin(dnu)/sqrt(mu*po);
	to   = go + sqrt(ao^3/mu)*(dEo-sin(dEo));

	a1   = m*k*p1/( (2*m-l*l)*p1^2 + 2*k*l*p1 - k^2);
	f1   = 1 - r2/p1 *(1-cos(dnu));
	dE1  = acos(1-r1/a1*(1-f1));
	g1   = r1*r2*sin(dnu)/sqrt(mu*p1);
	t1   = g1 + sqrt(a1^3/mu)*(dE1-sin(dE1));

	dt   = t1 - tx;
	while ( abs(dt) > 0.0000001 )
		p = p1 - dt*(p1-po)/(t1-to);
		po = p1; to = t1;
		p1 = p;
		a1   = m*k*p1/( (2*m-l*l)*p1^2 + 2*k*l*p1 - k^2);
		f1   = 1 - r2/p1 *(1-cos(dnu));
		dE1  = acos(1-r1/a1*(1-f1));
		g1   = r1*r2*sin(dnu)/sqrt(mu*p1);
		t1   = g1 + sqrt(a1^3/mu)*(dE1-sin(dE1));
		dt   = t1 - tx
	end;
	a = a1; dE = dE1;
