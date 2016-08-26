function [T] = ROT(ang,num)

if num == 1;
   T=[1, 0, 0; 0, cos(ang), sin(ang); 0, -sin(ang), cos(ang)];
elseif num == 2;
   T=[cos(ang), 0, -sin(ang); 0, 1, 0; sin(ang), 0, cos(ang)];
elseif num == 3;
   T=[cos(ang), sin(ang), 0; -sin(ang), cos(ang), 0; 0, 0, 1];
end
