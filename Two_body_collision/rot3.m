% Rotation about the Z axis
%
% rot3(a,gamma) where "a" is a cartesian vector and gamma is the
%               angle of rotation in radians

function [b] = rot3(a,gamma)


        b(1) = a(1)*cos(gamma) + a(2)*sin(gamma) + 0   ;
        b(2) =-a(1)*sin(gamma) + a(2)*cos(gamma) + 0   ;
        b(3) = 0               + 0               + a(3);
        b    = b';
