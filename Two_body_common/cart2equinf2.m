function [a,P1,P2,Q1,Q2,l]=cart2equinf2(x,y,z,xd,yd,zd)
%% this fuction converts cart (x,y,z,xd,yd,zd) to equinoctial (a,P1,P2,Q1,Q2,l) 
r0=[x;y;z];
v0=[xd;yd;zd];
mu=398601.2;
h=cross(r0,v0);
iz=[0,0,1]';
ix=[1,0,0]';
ih=h/norm(h);
a=1/(2/norm(r0)-norm(v0)^2/mu);
p=norm(h)^2/mu;
e=sqrt(1-p/a);
% i=acos(dot(ih,iz));
i=arc_angles('acos',dot(ih,iz));

ev=(1/mu*cross(v0,h)-1/norm(r0)*r0);
ie=(1/mu*cross(v0,h)-1/norm(r0)*r0)/norm(1/mu*cross(v0,h)-1/norm(r0)*r0);
ip=cross(ih,ie);
optss=optimset('Display','off','TolFun', 1.0e-15, 'TolX',1.0e-15);
in=fsolve(@(x)[dot([x(1);x(2);x(3)],ih);dot([x(1);x(2);x(3)],iz);1-norm([x(1);x(2);x(3)])],ie',optss);



Omg=arc_angles('acos',dot(ix,in));

omg=arc_angles('acos',dot(ie,in));
wb=Omg+omg;

% f0=acos(dot(ie,r0/norm(r0)));
f0=arc_angles('acos',dot(ie,r0/norm(r0)));

% E0=acos((e+cos(f0))/(1+e*cos(f0)));
E0=arc_angles('acos',(e+cos(f0))/(1+e*cos(f0)));
% if E0<0
%     E0=E0+2*pi;
% end
% if E0>0
%     E0=E0-2*pi;
% end
% M0=E0-e*sin(E0);
% b=sqrt(a^2*(1-e^2));
% P=2*pi*sqrt(a^3/mu);

l=wb+E0-e*sin(E0);
P1=e*sin(wb);
P2=e*cos(wb);
Q1=tan(i/2)*sin(Omg);
Q2=tan(i/2)*cos(Omg);

if isnan(P1)==1
    keyboard
end
