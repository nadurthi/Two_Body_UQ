% [a,phi] = R2aphi(R)
% for a given rotation matrix R,
% returns the Euler axis a (a 3x1 matrix)
% and Euler principal angle phi
% SEE ALSO aphi2R, aphi2q
function [a,phi] = R2aphi(R)
phi = acos((trace(R)-1)/2);
acr = (R'-R)/2/sin(phi);
a = [acr(3,2); acr(1,3); acr(2,1)];
a = a/norm(a);