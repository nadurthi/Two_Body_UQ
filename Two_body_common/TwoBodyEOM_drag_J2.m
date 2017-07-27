function rdot=TwoBodyEOM(t,r,Cd,msat,area);

%%
% This function returns the acceleration for a satellite orbiting in low
% Earth Orbit under J2 and Atmospheric Drag. 
%%

mu=398600.4415; %km^3/sec^2

Re=6378.1363; %km

rsat = r(1:3); rsat = rsat(:);
vsat = r(4:6); vsat = vsat(:);

ad_J2 = ComputeJ2(rsat,Re,mu); % J2 Calculation
adrag = ComputeAtmDarg(rsat,vsat,Re,Cd,area,msat); % Drag Calculation


R=norm(r(1:3));


rdot(1,1)=r(4);
rdot(2,1)=r(5);
rdot(3,1)=r(6);
rdot(4,1)=-mu*r(1)/R^3+ad_J2(1)+adrag(1);
rdot(5,1)=-mu*r(2)/R^3+ad_J2(2)+adrag(2);
rdot(6,1)=-mu*r(3)/R^3+ad_J2(3)+adrag(3);

function ad_J2 = ComputeJ2(rsat,Re,mu)

%%
% J2 Calculation
%%

J2=0.001082616; 

R = norm(rsat);

C11=rsat(1)/R;C12=rsat(2)/R;C13=rsat(3)/R;

adx=-1.5*J2*(mu/R^2)*(Re/R)^2*(1-5*C13^2)*C11;
ady=-1.5*J2*(mu/R^2)*(Re/R)^2*(1-5*C13^2)*C12;
adz=-1.5*J2*(mu/R^2)*(Re/R)^2*(3-5*C13^2)*C13;

ad_J2 = [adx;ady;adz];

function adrag = ComputeAtmDarg(rsat,vsat,Re,Cd,area,msat)


%%
% Atmospheric Drag Calculation
%%


Bal = Cd*area/msat; % Ballistic coefficient
R = norm(rsat);



omega_earth = 2*pi/(24*60*60); % Earth's rotation speed in radians per second
vatmos = [-omega_earth*rsat(2); omega_earth*rsat(1);0];   % km/sec
vrel = vsat-vatmos;  % relative velocity of satellite in km/sec
mvrel = norm(vrel);

rho = atmData(R-Re); %kg/m^3
rho = rho/(1e-3)^3; %kg/km^3

adrag =  -1/2*(Bal*rho*mvrel^2)*(vrel/mvrel);

