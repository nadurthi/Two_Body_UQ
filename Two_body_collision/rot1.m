% Rotation about the X axis
%
% rot1(a,alpha) where "a" is a cartesian vector and alpha is the
%               angle of rotation in radians

function [b] = rot1(a,alpha)


        b(1) = a(1) + 0               + 0              ;
        b(2) = 0    + a(2)*cos(alpha) + a(3)*sin(alpha);
        b(3) = 0    - a(2)*sin(alpha) + a(3)*cos(alpha);
        b    = b';
