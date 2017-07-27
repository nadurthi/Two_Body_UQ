function [x,y,z,xd,yd,zd]=equin2cart(a,P1,P2,Q1,Q2,l)
%% this fuction converts equinoctial (a,P1,P2,Q1,Q2,l) to cart (x,y,z,xd,yd,zd)
mu=398601.2;
wb=atan2(P1,P2);
if wb>2*pi
    wb=wb-2*pi;
end
if wb<0
    wb=wb+2*pi;
end

Omg=atan2(Q1,Q2);
if Omg>2*pi
    Omg=Omg-2*pi;
end
if Omg<0
    Omg=Omg+2*pi;
end

omg=wb-Omg;
if omg>2*pi
    omg=omg-2*pi;
end
if omg<0
    omg=omg+2*pi;
end

e=sqrt(P1^2+P2^2);
i1=2*atan2(sqrt(Q1^2+Q2^2),1);
i2=2*atan2(-sqrt(Q1^2+Q2^2),1);
if i1>2*pi
    i1=i1-2*pi;
end
if i1<0
    i1=i1+2*pi;
end
if i2>2*pi
    i2=i2-2*pi;
end
if i2<0
    i2=i2+2*pi;
end
if max(abs(Q1-tan(i1/2)*sin(Omg)),abs(Q2-tan(i1/2)*cos(Omg)))<=max(abs(Q1-tan(i2/2)*sin(Omg)),abs(Q2-tan(i2/2)*cos(Omg)))
    i=i1;
else
    i=i2;
end

optss=optimset('Display','off','TolFun', 1.0e-15, 'TolX',1.0e-15);
K=fsolve(@(K)K+P1*cos(K)-P2*sin(K)-l,0,optss);


% E=K-wb;
% M=E-e*sin(E);
b=a*sqrt(1-e^2);
p=a*(1-e^2);
h=sqrt(mu*p);
r=a*(1-P1*sin(K)-P2*cos(K));

abyab=a/(a+b);
sinL=(a/r)*((1-abyab*P2^2)*sin(K)+abyab*P1*P2*cos(K)-P1);
cosL=(a/r)*((1-abyab*P1^2)*cos(K)+abyab*P1*P2*sin(K)-P2);
L=atan2(sinL,cosL);
if L>2*pi
    L=L-2*pi;
end
if L<0
    L=L+2*pi;
end
f=L-wb;
% f=arc_angles('acos',(p-r)/(e*r));
% v=sqrt(mu*(2/r-1/a));
% keyboard
th=omg+f;
rv=r*[cos(Omg)*cos(th)-sin(Omg)*sin(th)*cos(i);...
    sin(Omg)*cos(th)+cos(Omg)*sin(th)*cos(i);...
    sin(th)*sin(i)]';

vv=-mu/h*[cos(Omg)*(sin(th)+e*sin(omg))+sin(Omg)*(cos(th)+e*cos(omg))*cos(i);...
          sin(Omg)*(sin(th)+e*sin(omg))-cos(Omg)*(cos(th)+e*cos(omg))*cos(i);...
          -(cos(th)+e*cos(omg))*sin(i)]';
    x=rv(1);
    y=rv(2);
    z=rv(3);
    xd=vv(1);
    yd=vv(2);
    zd=vv(3); 
