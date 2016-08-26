%    vs = ijk2sez(vi,L,th)
%    rotate from IJK to SEZ frame
%    vi = vector in inertial frame
%    L  = latitude
%    th = local sidereal time
%    vs = vector in SEZ frame
%    SEE ALSO sez2ijk
function vs = ijk2sez(vi,L,th)
	sL = sin(L);  cL = cos(L);
	sth = sin(th);  cth = cos(th);
	Rot =  [sL*cth -sth cL*cth;  sL*sth cth cL*sth;  -cL 0 sL];
	vs  =  Rot'*vi;
%
