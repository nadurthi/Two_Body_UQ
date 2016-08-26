function oe = oeofdt(oe,dt,mu)

%  oe = oeofdt(oe,dt,mu)
%
%  Author:  Capt Ralph Sandfry                   17 Feb 99
%
%  This function propogates the orbital elements in time.
%  This function is valid for non-circular, non-equatorial
%  orbits.
%
%  INPUTS:
%      oe      -  vector of the 6 orbital elements
%                 a (km),e,i,Om,om,nu (radians)
%      dt      -  time increment                          sec
%      mu      -  Gravitational Parameter                 km^3/s^2
%
%  OUTPUTS:
%      oe      -  updated orbital element vector
%
%  SEE ALSO rv2oe, oe2rv, tle2oe

twopi = 2*pi;

%  LOCALS:
%      I       -  counter
%      ai      -  initial semi-major axis                 km
%      af      -  final semi-major axis                   km
%      ei      -  initial eccentricity
%      ef      -  final eccentricity
%      i_i     -  initial inclination                     rad
%      i_f     -  final inclination                       rad
%      Om_i    -  initial longitude of ascending node     rad
%      Om_f    -  final longitude of ascending node       rad
%      om_i    -  initial argument of periapsis           rad
%      om_f    -  final argument of periapsis             rad
%      nu_i    -  initial true anomaly                    rad
%      nu_f    -  final true anomaly                      rad
%      Ei      -  initial eccentric anomaly               rad
%      Ef      -  final eccentric anomaly                 rad
%      Mi      -  initial mean anomaly                    rad
%      Mf      -  final mean anomaly                      rad
%      n       -  mean motion                             rad/s
%
%  CONSTANTS:
%      twopi   -  2*pi
%
%  COUPLING:
%      none
%
%  REFERENCES:
%      Understanding Space, Ch 8
%


for I=3:6                %  check for negative angles
   if oe(I)<0            %  and make 0 < angle < 2*pi
      oe(I)=oe(I)+twopi;
   end
end

ai = oe(1);
ei = oe(2);
i_i = oe(3);
Omi = oe(4);
omi = oe(5);
nui = oe(6);

%  Update first five orbital elements

af = ai;
ef = ei;
i_f = i_i;
omf = omi;
Omf = Omi;

% Calculate initial eccentric anomaly

Ei = acos((ei+cos(nui))/(1+ei*cos(nui)));

if nui > pi                    % half-plane check
   Ei = 2*pi - Ei;
end

Mi = Ei - ei*sin(Ei);          % initial mean anomaly

n = sqrt(mu/ai^3);             %  mean motion

Mf = Mi + n*dt;                
Mf = mod(Mf,twopi);            %  future mean anomaly

Ef = EofMe(Mf,ef);             %  future eccentric anomaly

nuf = acos((cos(Ef)-ef)/(1-ef*cos(Ef)));
if Ef > pi                     
   nuf = twopi - nuf;           %  future true anomaly with
end                             %  half plane check 
   
oe = [af;ef;i_f;Omf;omf;nuf];

