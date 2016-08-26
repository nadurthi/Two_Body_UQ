% q = R2q(R)
%     R is a rotation matrix
%     q is a quaternion (4 x 1) with q(1:3) the vector part
function q = R2q(R)
q4 = 0.5*sqrt(1+trace(R));
q13 = [ R(2,3)-R(3,2); R(3,1)-R(1,3); R(1,2)-R(2,1)]/4/q4;
q = [q13;q4];
q = q/norm(q);   % make sure it's unit length
