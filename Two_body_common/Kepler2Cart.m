function [r_v,v_v] = Kepler2Cart(a,e,i,w,O,f)
%	[r_v,v_v] = Kepler2Cart(a,e,ini_ang,Arg_Par,RAAN,True_Anomaly)
%  This function takes the 6 kepler orbit elements (all angles in radians)
%    and converts them to the radius and velocity vectors in the ECI
%    frame.  This code uses the equations given in Design 4 of the ASEN
%    2004 class (Also posted on this class's web site).  You must make the 
%    variable mu (G*Me) global in your main program to get this to work.  
global mu

%Compute the radius mag
r = a*(1-e^2)/(1+e*cos(f));

%Compute the specific angular momentum (usually h)
l = sqrt(mu*a*(1-e^2));

%Find the X-Y-Z cosines for r vector
al = w+f;%Arg of Latitude
x = cos(O)*cos(al)-sin(O)*sin(al)*cos(i);
y = sin(O)*cos(al)+cos(O)*sin(al)*cos(i);
z = sin(i)*sin(al);

%Calculate r vector
r_v = r*[x y z];

%Orbit parameter p (Not period)
p = a*(1-e^2);

%Find the velocity Cosines
d_v = [cos(O)*sin(al)+sin(O)*cos(al)*cos(i) ...
sin(O)*sin(al)-cos(O)*cos(al)*cos(i) -sin(i)*cos(al)];

%Find the velocity
v_v = r_v*l*e*sin(f)/r/p-l/r*d_v;
