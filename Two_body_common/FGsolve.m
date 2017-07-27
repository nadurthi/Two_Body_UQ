function  [T,Y,F,G,Ehatt] = FGsolve(x0,t)
% F and G series solution - analytic solution to the 2 Body Problem
% 
mue = 398601.2;  % Gravitational parameter of earth
tol = 1e-15;
R0 = x0(1:3,:); V0 = x0(4:6,:);
r0 = sqrt(R0'*R0); v0 = sqrt(V0'*V0);
a = 1/(2/r0 - v0^2/mue);
sig0 = R0'*V0/sqrt(mue);

Y(1,1:6)=x0'; T = t;
Ehat = 0;

for k = 2:length(t)
    dt = t(k)-t(1);
    Ehat0 = Ehat;
    [Ehat]=keplermod(dt, r0, a, sig0, Ehat0, tol);
    Ehatt(k,1)=Ehat;
    
    r = a + (r0 - a)*cos(Ehat) + sqrt(a)*sig0*sin(Ehat);
    % F and G functions
    F(k,1) = 1-a*(1-cos(Ehat))/r0;
    G(k,1) = dt + sqrt(a^3/mue)*(sin(Ehat)-Ehat);
    % Ft and Gt functions
    F(k,2) = -sqrt(mue*a)*sin(Ehat)/(r*r0);
    G(k,2) = 1+a*(cos(Ehat)-1)/r;
    
    Rt = F(k,1)*R0 + G(k,1)*V0;
    Vt = F(k,2)*R0 + G(k,2)*V0;
    Y(k,1:3) = Rt'; 
    Y(k,4:6) = Vt';
    
end

