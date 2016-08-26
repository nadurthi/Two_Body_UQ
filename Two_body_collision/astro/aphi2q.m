% q = aphi2q(a,phi)
%     a is the Euler axis (a 3 x 1 matrix)
%     phi is the Euler principal angle, a scalar
%     returns the quaternion
% SEE ALSO q2aphi
function q = aphi2q(a,phi)
[m,n]=size(a);
if (m==1 & n==3)
    a=a';
end
if ((m~=1 & n~=1) | (m~=3 & n~=3))
   disp('a is not a 3 x 1 matrix');
   q=NaN*ones(4,1);
   return;
end
c = cos(phi/2); s = sin(phi/2);
q = [a*s; c];