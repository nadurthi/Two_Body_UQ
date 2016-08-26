% R = q2R(q)
%     q = quaternion (4x1)
%     R = rotation matrix
% SEE ALSO R2q
function R = q2R(q)
R = (q(4)^2-q(1:3)'*q(1:3))*eye(3)+2*q(1:3)*q(1:3)'-2*q(4)*cr(q(1:3));