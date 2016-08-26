%  vi = sez2ijk(vs,L,th)  
%       rotate from SEZ to IJK frame
%  vs = vector in SEZ frame
%  L  = latitude
%  th = local sidereal time
%  vi = vector in IJK frame
function vi = sez2ijk(vs,L,th)
	sL = sin(L);  cL = cos(L);
	sth = sin(th);  cth = cos(th);
	Rot =  [sL*cth -sth cL*cth;  sL*sth cth cL*sth;  -cL 0 sL];
	vi  =  Rot*vs;