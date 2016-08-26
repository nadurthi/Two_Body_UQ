%   H = lineft(r1,r2,thickness,color)
%   draws a line of thickness and color from r1 to r2
%   r1,r2 = 3x1 position vectors
%   H = handle to line object
function H = lineft(r1,r2,thickness,color)
    if nargin==2, thickness=1; color='b';  end;
    if nargin==3, color='b'; end;
    X = [r1(1) r2(1)]';
    Y = [r1(2) r2(2)]';
    Z = [r1(3) r2(3)]';
    H = line(X,Y,Z);
    set(H,'LineWidth',thickness,'color',color);
