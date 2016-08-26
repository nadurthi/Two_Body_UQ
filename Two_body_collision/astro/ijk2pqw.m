%
%  ijk2pqw.m  r,v in ijk to pqw
%
%  [r,v] = ijk2pqw(r,v,mu)
%			r,v  expressed in  IJK  frame
%        returns r,v in PQW frame
%
function [rv,vv] = ijk2pqw(ri,vi,mu)
   oe=rv2oe(ri,vi,mu);
	a=oe(1); e=oe(2); i=oe(3); Om=oe(4); om=oe(5); nu=oe(6);
	p = a*(1-e*e);
	r = p/(1+e*cos(nu));
	rv = [r*cos(nu); r*sin(nu); 0];			% in PQW frame
	vv = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0];
