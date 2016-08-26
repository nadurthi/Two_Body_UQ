%  razel2rv.m  RhoAzEl to  r,v  conversion
%
%  [r,v] = razel2rv(rhoazel,th,L,R,w)
%			rhoazel = [rho rhod az azd el eld]
%                                  km, km/s, rad, rad/s
%			th = local sidereal time (rad)
%			L  = site latitude (rad)
%			R  = radius of earth (km)
%			w  = angular velocity of earth (rad/s)
%        r  = position vector
%        v = velocity vector
function [ri,vi] = razel2rv(rhoazel,th,L,R,w)
	rho =  rhoazel(1);   rhod = rhoazel(2);
	az  =  rhoazel(3);   azd  = rhoazel(4);
	el  =  rhoazel(5);   eld  = rhoazel(6);
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
	rhov=  [rs; re; rz];
%
%  components of rho^circle in SEZ frame
%
	rsd = -rhod*ce*ca+rho*se*eld*ca+rho*ce*sa*azd;
	red =  rhod*ce*sa-rho*se*eld*sa+rho*ce*ca*azd;
	rzd =  rhod*se+rho*ce*eld;
	rhoc=  [rsd; red; rzd];
%
%  position of site in SEZ frame
%
	Rv  =  R*[0;0;1];
%
%  position of vehicle in SEZ frame
%
	rv  =  Rv + rhov;
%
%  rotate into IJK frame
%
	Rot =  [sL*cth -sth cL*cth;  sL*sth cth cL*sth;  -cL 0 sL];
	ri  =  Rot*rv;
%
%  velocity in IJK frame
%
	vi  =  Rot*rhoc + cross([0;0;w],ri);
