clear all; close all; clc;
global mue re;

re = 6372.797;   % Mean earth radius (Km)
mue = 398601.2;  % Gravitational parameter of earth (Km^3/Sec^2)


r0 = [7000; 0; 0];
rdot0 = [0; -1.03741; 7.977129];

[a0, e0, E0, w0, i0, Om0] = XYZ2OE(r0,rdot0,mue)


M0 = E0 - e0*sin(E0);

Nsat = 20;

% uniform distribution parameters
M0a = 0; M0b = pi/5;

M0sat = M0 + M0a + (M0b-M0a)*rand(Nsat,1);

for j = 1:Nsat
    E0sat = kepler(M0sat(j), e0, 0, 1e-12);
    [X1, Xdot1] = OE2XYZ(a0, e0, E0sat, w0, i0, Om0, mue);
    Xhist(:,j) = X1;
    Xdothist(:,j) = Xdot1;
end


figure
plot3(Xhist(1,:),Xhist(2,:),Xhist(3,:),'*')

