%
%  oep2rv.m  Orbital Elements to r,v for parabolic orbit
%
%  [r,v] = oep2rv(oep,mu)
%			oep = [p e i Om om nu]
%			r,v  expressed in  IJK  frame
%
function [ri,vi] = oep2rv(oep,mu)
	p=oep(1); e=oep(2); i=oep(3); Om=oep(4); om=oep(5); nu=oep(6);
	r = p/(1+e*cos(nu));
	rv = [r*cos(nu); r*sin(nu); 0];			% in PQW frame
	vv = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0];
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
