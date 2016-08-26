%
%  Rijk2pqw.m  Orbital Elements to Rpi
%
%  Rpi = Rijk2pqw(oe)
%			oe = [a e i Om om nu]
%
function Rpi = Rijk2pqw(oe)
	i=oe(3); Om=oe(4); om=oe(5); 
	cO = cos(Om);  sO = sin(Om);
	co = cos(om);  so = sin(om);
	ci = cos(i);   si = sin(i);
	Rip  = [cO*co-sO*so*ci  -cO*so-sO*co*ci  sO*si;
		  sO*co+cO*so*ci  -sO*so+cO*co*ci -cO*si;
        so*si            co*si           ci];
   Rpi=Rip';
