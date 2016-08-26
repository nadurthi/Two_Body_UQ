%  function R = Roi(Om,inc,u)
%  returns the rotation matrix from inertial to
%  orbital frame, defined by
%    Om  = right ascension of the ascending node
%    inc = inclination
%    u   = argument of latitude
function R = Roi(Om,inc,u)
sO = sin(Om);  cO = cos(Om);
si = sin(inc); ci = cos(inc);
su = sin(u);   cu = cos(u);
R = [-su*cO-cu*ci*sO -su*sO+cu*ci*cO cu*si;
     -si*sO si*cO -ci;
     -cu*cO+su*ci*sO -cu*sO-su*ci*cO -su*si];
   