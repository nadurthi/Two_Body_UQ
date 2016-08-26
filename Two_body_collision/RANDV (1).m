function [posijk,velijk] = randv(a,e,i,Omega,w,nu)

mu = 398600.4;	%km^3/s^2

p = a*(1-e^2);

i=i*pi/180;
Omega=Omega*pi/180;
w=w*pi/180;
nu=nu*pi/180;

Rpqw = [p*cos(nu)/(1+e*cos(nu)); p*sin(nu)/(1+e*cos(nu)); 0];

Vpqw = [-sqrt(mu/p)*sin(nu); sqrt(mu/p)*(e+cos(nu)); 0];

posijk = ROT(-Omega,3)*ROT(-i,1)*ROT(-w,3)*Rpqw;
velijk = ROT(-Omega,3)*ROT(-i,1)*ROT(-w,3)*Vpqw;

