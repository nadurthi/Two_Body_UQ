%  razel2r.m  RhoAzEl to  r  conversion
%
%  r = razel2r(rhoazel,th,L,R)
%			rhoazel = [rho az el]
%			th = local sidereal time
%			L  = site latitude
%			R  = radius of earth
%        r  = position vector in IJK frame
%  SEE ALSO razel2rv, rv2oe, oe2rv
function ri = razel2r(rhoazel,th,L,R)
	rho =  rhoazel(1);
	az  =  rhoazel(2);
	el  =  rhoazel(3);
	ce  =  cos(el);      se   = sin(el);
	ca  =  cos(az);      sa   = sin(az);
	cL  =  cos(L);       sL   = sin(L);
	cth =  cos(th);      sth  = sin(th);
%
%  componenets of rho in SEZ frame
%
	rs  = -rho*ce*ca;
	re  =  rho*ce*sa;
	rz  =  rho*se;
	rhov=  [rs; re; rz]
%
%  position of site in SEZ frame
%
	Rv  =  R*[0;0;1];
%
%  position of vehicle in SEZ frame
%
	rv  =  Rv + rhov
%
%  rotate into IJK frame
%
	Rot =  [sL*cth -sth cL*cth;  sL*sth cth cL*sth;  -cL 0 sL]
	ri  =  Rot*rv;
