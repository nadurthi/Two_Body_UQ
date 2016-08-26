function deqOrb=diff_equinoctial(t,eqOrb)
a=eqOrb(1);
P1=eqOrb(2);
P2=eqOrb(3);
Q1=eqOrb(4);
Q2=eqOrb(5);
l=eqOrb(6);


Cd = 3; % Assuming satellite to be sphere
msat = 1000; %in kg
area = 100*(1e-3)^2; % in km^2
mu=398601.2;
Re=6378.165; %km

[x,y,z,xd,yd,zd]=equin2cart(a,P1,P2,Q1,Q2,l);
rv=[x;y;z];
vv=[xd;yd;zd];

% e=sqrt(P1^2+P2^2);
% wb=atan2(P1,P2);
b=a*sqrt(1-P1^2-P2^2);

optss=optimset('Display','off');
K=fsolve(@(K)K+P1*cos(K)-P2*sin(K)-l,0,optss);
r=a*(1-P1*sin(K)-P2*cos(K));

abyab=a/(a+b);
sinL=(a/r)*((1-abyab*P2^2)*sin(K)+abyab*P1*P2*cos(K)-P1);
cosL=(a/r)*((1-abyab*P1^2)*cos(K)+abyab*P1*P2*sin(K)-P2);
L=atan2(sinL,cosL);

n=sqrt(mu/a^3);
h=n*a*b;
pbyr=1+P1*sinL+P2*cosL;
rbyh=h/(mu*(1+P1*sinL+P2*cosL));

ad_J2 = ComputeJ2(rv,Re,mu); % J2 Calculation
adrag = ComputeAtmDarg(rv,vv,Re,Cd,area,msat); % Drag Calculation

[a,e,i,omg,Omg,M]=equin2orbital(a,P1,P2,Q1,Q2,l);

l1=cos(Omg)*cos(omg)-sin(Omg)*sin(omg)*cos(i);
l2=-cos(Omg)*sin(omg)-sin(Omg)*cos(omg)*cos(i);
l3=sin(Omg)*sin(i);
m1=sin(Omg)*cos(omg)+cos(Omg)*sin(omg)*cos(i);
m2=-sin(Omg)*sin(omg)+cos(Omg)*cos(omg)*cos(i);
m3=-cos(Omg)*sin(i);
n1=sin(omg)*sin(i);
n2=cos(omg)*sin(i);
n3=cos(i);
R=[l1,l2,l3;m1,m2,m3;n1,n2,n3];
wb=Omg+omg;
E=K-wb;

Rf=[(a/r)*(cos(E)-e),-(b/r)*sin(E),0;(b/r)*sin(E),(a/r)*(cos(E)-e),0;0,0,1];
aa=inv(R*Rf)*(ad_J2+adrag);

adr=1*aa(1);
adth=1*aa(2);
adh=1*aa(3);

dadt=(2*a^2/h)*((P2*sinL-P1*cosL)*adr+pbyr*adth);

dP1dt=rbyh*(-pbyr*cosL*adr+(P1+(1+pbyr)*sinL)*adth-P2*(Q1*cosL-Q2*sinL)*adh);

dP2dt=rbyh*(pbyr*sinL*adr+(P2+(1+pbyr)*cosL)*adth+P1*(Q1*cosL-Q2*sinL)*adh);

dQ1dt=0.5*rbyh*(1+Q1^2+Q2^2)*sinL*adh;

dQ2dt=0.5*rbyh*(1+Q1^2+Q2^2)*cosL*adh;

dldt=n-rbyh*((abyab*pbyr*(P1*sinL+P2*cosL)+2*b/a)*adr+abyab*(1+pbyr)*(P1*cosL-P2*sinL)*adth+(Q1*cosL-Q2*sinL)*adh);

deqOrb=[dadt;dP1dt;dP2dt;dQ1dt;dQ2dt;dldt];
% keyboard
t
end

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

end

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
% return;
end