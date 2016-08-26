% xc = cr(x)
%      x = 3x1 matrix
%      xc = skew symmetric 3x3 matrix
%           [ 0   -x(3),   x(2)
%             x(3)  0,    -x(1)
%            -x(2)  x(1)    0  ]
function xc = cr(x)
xc = [ 0   -x(3),   x(2)
       x(3)  0,    -x(1)
      -x(2)  x(1)    0 ];


