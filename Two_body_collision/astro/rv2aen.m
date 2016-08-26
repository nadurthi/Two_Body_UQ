%
%  rv2aen.m  r,v,mu -> a, e, and nu
%
%  aen = rv2aen(rv,vv,mu)
%        rv  = [X Y Z]'    vv = [v1 v2 v3]'   in IJK frame
%        aen = [a e nu0]
%
function aen = rv2aen(rv,vv,mu)
%
% first do some preliminary calculations
%
	hv = cr(rv)*vv;
	h2 = (hv'*hv);
	v2 = (vv'*vv);
	r  = sqrt(rv'*rv);
	ev = 1/mu *( (v2-mu/r)*rv - (rv'*vv)*vv );
	p  = h2/mu;
%
% now compute the required oe's
%
	e  = sqrt(ev'*ev);		% eccentricity
	a  = p/(1-e*e);			% semimajor axis
	
	nu = acos(ev'*rv/e/r);		% true anomaly
	if ( rv'*vv < 0 )		% fix quadrant
		nu = 2*pi-nu;
	end;
	aen = [a e nu];			% assemble "vector"
