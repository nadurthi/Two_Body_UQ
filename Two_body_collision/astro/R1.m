% R = R1(angle)
% returns elementary "1" rotation matrix
function Rotmat = R1(angle)
   c = cos(angle); s = sin(angle);
   Rotmat = [ 1 0 0; 0 c s; 0 -s c ];
