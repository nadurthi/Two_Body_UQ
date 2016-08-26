%  rv2oe.m  r,v,mu -> orbital elements
%
%  oe = rv2oe(rv,vv,mu)
%       rv  = [X Y Z]'    vv = [v1 v2 v3]'   in IJK frame
%       oe = [a e i Omega omega nu0]
%
%  SEE ALSO oe2rv, tle2oe
function oe = rv2oe(rv,vv,mu)
%
% first do some preliminary calculations
%
	K  = [0;0;1];
	hv = cross(rv,vv);
	nv = cross(K,hv);
	n  = sqrt(nv'*nv);
	h2 = (hv'*hv);
	v2 = (vv'*vv);
	r  = sqrt(rv'*rv);
   ev = 1/mu *( (v2-mu/r)*rv - (rv'*vv)*vv );
	p  = h2/mu;
%
% now compute the oe's
%
	e  = sqrt(ev'*ev);		% eccentricity
	a  = p/(1-e*e);			% semimajor axis
	i  = acos(hv(3)/sqrt(h2));	% inclination
   
   argacos=nv(1)/n;
   if (abs(argacos)>1)
      argacos=sign(argacos);
   end
   Om = acos(argacos);		% RAAN
	if ( nv(2) < 0 )		% fix quadrant
		Om = 2*pi-Om;
   end;
   
   argacos=nv'*ev/n/e;
	if (abs(argacos)>1)
      argacos=sign(argacos);
   end
   om = acos(argacos);		% arg of periapsis
	if ( ev(3) < 0 )		% fix quadrant
		om = 2*pi-om;
   end;
   
   argacos = ev'*rv/e/r;
   if (abs(argacos)>1)
      argacos=sign(argacos);
   end
   nu = acos(argacos);		% true anomaly
	if ( rv'*vv < 0 )		% fix quadrant
		nu = 2*pi-nu;
	end;
	oe = [a e i Om om nu];		% assemble "vector"
