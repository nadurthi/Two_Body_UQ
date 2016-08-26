function Y = equinoxToCart(X)
a = X(1);
h = X(2);
k = X(3);
p = X(4);
q = X(5);
l = X(6);
MU = 398601.2;%398600.4418;
syms x 

K = solve('K + h*cos(K)-k*sin(K) = l');
r1 = a*(1 - h*sin(K) - k*cos(K));
e = sqrt(h^2+k^2);
b = a*sqrt(1 - e^2);
x1 = double(a/r1*((1-a/(a+b)*k^2)*sin(K) + a/(a+b)*h*k*cos(K)-h));
y1 = double(a/r1*((1-a/(a+b)*h^2)*cos(K) + a/(a+b)*h*k*sin(K)-k));
L = atan2(x1,y1);
r = a*(1-h^2-k^2)/(1+k*cos(L)+h*cos(L));
s = sqrt(1+p^2+q^2);
w = 1 + k*cos(L) + h*sin(L);
p1 = a*(1-e^2);
alpha = real(sqrt(p^2 - q^2));

cartX = r/s^2*(cos(L) + alpha^2*cos(L) + 2*p*q*sin(L));
cartY = r/s^2*(sin(L) + alpha^2*sin(L) + 2*p*q*cos(L));
cartZ = 2*r/s^2*(q*sin(L) - p*cos(L));
cartXDot = -1/s^2*sqrt(MU/p1)*(sin(L) + alpha^2*sin(L) - 2*p*q*cos(L) + h - 2*k*p*q + alpha^2*h);
cartYDot = -1/s^2*sqrt(MU/p1)*(-cos(L) + alpha^2*cos(L) + 2*p*q*sin(L) - k + 2*h*p*q + alpha^2*k);
cartZDot = 2/s^2*sqrt(MU/p1)*(q*cos(L) + p*sin(L) + p*k + h*p);

Y = [double(cartX);double(cartY);double(cartZ);double(cartXDot);double(cartYDot);double(cartZDot)]';