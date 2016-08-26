%  oee2rv.m  Orbital Elements (E) to r,v
%
%  [r,v] = oee2rv(oe,mu)
%			oe = [a e i Om om E]
%			r,v  expressed in  IJK  frame
%
%  SEE ALSO oe2rv, rv2oe, oee2r, tle2oe
function [ri,vi] = oee2rv(oe,mu)
	a=oe(1); e=oe(2); i=oe(3); Om=oe(4); om=oe(5); E=oe(6);
   cE = cos(E);   sE = sin(E);
	p = a*(1-e*e);
	r = a*(1-e*cE);
        cnu = (e-cE)/(e*cE-1);
        snu = a*sqrt(1-e*e)/r*sE;
	rv = [r*cnu; r*snu; 0];			% in PQW frame
	vv = sqrt(mu/p)*[-snu; e+cnu; 0];
%
%	now rotate
%
	cO = cos(Om);  sO = sin(Om);
	co = cos(om);  so = sin(om);
	ci = cos(i);   si = sin(i);
	R  = [cO*co-sO*so*ci  -cO*so-sO*co*ci  sO*si;
		  sO*co+cO*so*ci  -sO*so+cO*co*ci -cO*si;
		  so*si            co*si           ci];
	ri = R*rv;
	vi = R*vv;
