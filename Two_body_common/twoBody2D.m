function rdot = twoBody2D(t,x)
r =x;
Cd = 3; % Assuming satellite to be sphere
msat = 1000; %in kg
area = 100*(1e-3)^2; % in km^2
mu=398601.2; %km^3/sec^2

Re=6378.165; %km

rsat = r(1:2); rsat = rsat(:);
vsat = r(3:4); vsat = vsat(:);

% ad_J2 = 1*ComputeJ2(rsat,Re,mu); % J2 Calculation
% adrag = 1*ComputeAtmDarg(rsat,vsat,Re,Cd,area,msat); % Drag Calculation
% ad_J2=[0,0,0];
% adrag=[0,0,0];

R=norm(r(1:2));

rdot(1,1)=r(3);
rdot(2,1)=r(4);
rdot(3,1)=-mu*r(1)/R^3;
rdot(4,1)=-mu*r(2)/R^3;



end
