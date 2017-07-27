clear all; close all; clc;
global mue re;

re = 6372.797;   % Mean earth radius (Km)
mue = 398601.2;  % Gravitational parameter of earth (Km^3/Sec^2)


r0 = [7000; 0; 0];
rdot0 = [0; -1.03741; 7.977129];

[a0, e0, E0, w0, i0, Om0] = XYZ2OE(r0,rdot0,mue)


M0 = E0 - e0*sin(E0);

Nsat = 50;
n1 = 11;
% generation of initial conditions 
Xhist(1:6,Nsat) = 0;
for j = 1:Nsat 
    if(j < 47)
        OEset(:,j) = [a0; e0; E0; w0; i0; Om0+(j-1)*2*pi/Nsat];
        Om0+j*pi/(Nsat)
        [Xhist(1:3,j), Xhist(4:6,j)]= OE2XYZ(OEset(1,j), OEset(2,j), OEset(3,j), ...
                            OEset(4,j), OEset(5,j), OEset(6,j), mue);
    else
         
         switch j
             case 48 
             e0 = 0.15;
             OEset(:,j) = [a0; e0; E0; w0; i0; Om0+(j-1)*2*pi/Nsat];
             [Xhist(1:3,j), Xhist(4:6,j)]= OE2XYZ(OEset(1,j), OEset(2,j), OEset(3,j), ...
                            OEset(4,j), OEset(5,j), OEset(6,j), mue);
             case 49
             e0 = 0.18;
             OEset(:,j) = [a0; e0; E0; w0; i0; Om0+(j-1)*2*pi/Nsat];
             [Xhist(1:3,j), Xhist(4:6,j)]= OE2XYZ(OEset(1,j), OEset(2,j), OEset(3,j), ...
                            OEset(4,j), OEset(5,j), OEset(6,j), mue);    
             case 50
             e0 = 0.2;
             OEset(:,j) = [a0; e0; E0; w0; i0; Om0+(j-1)*2*pi/Nsat];
             [Xhist(1:3,j), Xhist(4:6,j)]= OE2XYZ(OEset(1,j), OEset(2,j), OEset(3,j), ...
                            OEset(4,j), OEset(5,j), OEset(6,j), mue);
         end
    end
end

% simulations
% orbital period
tp = 2*pi*sqrt(a0^3/mue);
t = linspace(0,2*tp,1000); tspan = t;

for j = 1:Nsat
    xb0 = [Xhist(1:3,j);Xhist(4:6,j)];
    % F&G series solution
    [T1,Y1,F,G,Ehatt] = FGsolve(xb0,t);
    
    xh(:,j) = Y1(:,1); yh(:,j) = Y1(:,2);
    zh(:,j) = Y1(:,3); xdoth(:,j) = Y1(:,4);
    ydoth(:,j) = Y1(:,5); zdoth(:,j) = Y1(:,6);
end


%% plot all the orbits. 

[Xe,Ye,Ze] = sphere(20);

figure
plot3(xh(:,1), yh(:,1), zh(:,1)); hold on;
    for j = 2: Nsat
        plot3(xh(:,j),yh(:,j),zh(:,j))
    end
mesh(Xe*re,Ye*re,Ze*re)
axis equal
