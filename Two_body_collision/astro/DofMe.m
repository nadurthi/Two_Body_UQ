%  function D = DofMe(M,p,tol)
%
%  this function computes D as a function of M and p
%
function D = DofMe(M,p,tol)
	Dn  = M
   Dn1 = Dn - (p*Dn+Dn^3/3-M)/(p+Dn^2)
   while ( abs(Dn1-Dn) > tol )
		Dn = Dn1;
		   Dn1 = Dn - (p*Dn+Dn^3/3-M)/(p+Dn^2)
	end;
	D = Dn1;
