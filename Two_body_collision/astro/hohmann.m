% [dv1,dv2,dt]=hohmann(r1,r2,mu)
% computes the two delta v's for a hohmann transfer from r1 to r2 
% and computes the time of flight
% for a planet with gravitational parameter mu
function [dv1,dv2,dt]=hohmann(r1,r2,mu)
vc1=sqrt(mu/r1);
vc2=sqrt(mu/r2);
at=(r1+r2)/2;
Et=-mu/2/at;
v1=sqrt(2*(Et+mu/r1));
v2=sqrt(2*(Et+mu/r2));
dv1=v1-vc1;
dv2=vc2-v2;
dt=pi*sqrt(at^3/mu);


