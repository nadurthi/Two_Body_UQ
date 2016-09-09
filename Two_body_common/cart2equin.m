function [a,P1,P2,Q1,Q2,l]=cart2equin(x,y,z,xd,yd,zd)
%% this fuction converts cart (x,y,z,xd,yd,zd) to equinoctial (a,P1,P2,Q1,Q2,l) 
r0=[x;y;z];
v0=[xd;yd;zd];
mu=398601.2;
hv=cross(r0,v0);
iz=[0,0,1]';
ix=[1,0,0]';
iy=[0,1,0]';
ih=hv/norm(hv);
a=1/(2/norm(r0)-norm(v0)^2/mu);
p=norm(hv)^2/mu;
e=sqrt(1-p/a);
ev=1/mu*(cross(v0,hv)-(mu/norm(r0))*r0);
ie=ev/norm(ev);
ip=cross(ih,ie);
h=norm(hv);
[Omg,omg,i]=solve_rot_mat(ix,iy,iz,ie,ip,ih);

% i=acos(dot(ih,iz));
% i=arc_angles('acos',dot(ih,iz));

% ev=(1/mu*cross(v0,h)-1/norm(r0)*r0);
% ie=(1/mu*cross(v0,h)-1/norm(r0)*r0)/norm(1/mu*cross(v0,h)-1/norm(r0)*r0);
% ip=cross(ih,ie);
% optss=optimset('Display','off','TolFun', 1.0e-15, 'TolX',1.0e-15);
% % in=fsolve(@(x)[dot([x(1);x(2);x(3)],ih);dot([x(1);x(2);x(3)],iz);1-norm([x(1);x(2);x(3)])],ie',optss);
% z1=iz(1);
% z2=iz(2);
% z3=iz(3);
% h1=ih(1);
% h2=ih(2);
% h3=ih(3);
% % 
% if h2*z1-h1*z2~=0
% n3=-abs(h2*z1-h1*z2)/sqrt((h3^2*(z1^2+z2^2)-2*h1*h3*z1*z3-2*h2*z2*(h1*z1+h3*z3)+h2^2*(z1^2+z3^2)+h1^2*(z2^2+z3^2)));
% n1=n3*(h3*z2-h2*z3)/(h2*z1-h1*z2);
% n2=n3*(h3*z1-h1*z3)/(-h2*z1+h1*z2);
% in1=[n1;n2;n3];
% 
% n3=abs(h2*z1-h1*z2)/sqrt((h3^2*(z1^2+z2^2)-2*h1*h3*z1*z3-2*h2*z2*(h1*z1+h3*z3)+h2^2*(z1^2+z3^2)+h1^2*(z2^2+z3^2)));
% n1=n3*(h3*z2-h2*z3)/(h2*z1-h1*z2);
% n2=n3*(h3*z1-h1*z3)/(-h2*z1+h1*z2);
% in2=[n1;n2;n3];
% if norm(in1-ie)>norm(in2-ie)
%     in=in2;
% else
%     in=in1;
% end
% 
% elseif h3*z1-h1*z3~=0
% n2=-1/sqrt(1+((h2*z1-h1*z2)^2/(h3*z1-h1*z3)^2)+((h3*z2-h2*z3)^2/(h3*z1-h1*z3)^2));
% n1=n2*(h3*z2-h2*z3)/(-h3*z1+h1*z3);
% n3=n2*(h2*z1-h1*z2)/(-h3*z1+h1*z3);
% in1=[n1;n2;n3];
% 
% n2=1/sqrt(1+((h2*z1-h1*z2)^2/(h3*z1-h1*z3)^2)+((h3*z2-h2*z3)^2/(h3*z1-h1*z3)^2));
% n1=n2*(h3*z2-h2*z3)/(-h3*z1+h1*z3);
% n3=n2*(h2*z1-h1*z2)/(-h3*z1+h1*z3);
% in2=[n1;n2;n3];
% if norm(in1-ie)>norm(in2-ie)
%     in=in2;
% else
%     in=in1;
% end
% 
% elseif h3*z2-h2*z3~=0
% n1=-abs(h3*z2-h2*z3)/sqrt(h3^2*(z1^2 + z2^2) - 2 *h1* h3 *z1* z3 - 2* h2* z2* (h1* z1 + h3* z3) +  h2^2 *(z1^2 + z3^2) + h1^2* (z2^2 + z3^2));
% n2=n1*(h3*z1-h1*z3)/(-h3*z2+h2*z3);
% n3=n1*(h2*z1-h1*z2)/(h3*z2-h2*z3);
% in1=[n1;n2;n3];
% 
% n1=abs(h3*z2-h2*z3)/sqrt(h3^2*(z1^2 + z2^2) - 2 *h1* h3 *z1* z3 - 2* h2* z2* (h1* z1 + h3* z3) +  h2^2 *(z1^2 + z3^2) + h1^2* (z2^2 + z3^2));
% n2=n1*(h3*z1-h1*z3)/(-h3*z2+h2*z3);
% n3=n1*(h2*z1-h1*z2)/(h3*z2-h2*z3);
% in2=[n1;n2;n3];
% if norm(in1-ie)>norm(in2-ie)
%     in=in2;
% else
%     in=in1;
% end
% else
%     in=fsolve(@(x)[dot([x(1);x(2);x(3)],ih);dot([x(1);x(2);x(3)],iz);1-norm([x(1);x(2);x(3)])],ie',optss);
% end

% keyboard
% Omg=arc_angles('acos',dot(ix,in));

% omg=arc_angles('acos',dot(ie,in));
wb=Omg+omg;

% f0=acos(dot(ie,r0/norm(r0)));
% f0=arc_angles('acos',dot(ie,r0/norm(r0)));
ir=r0/norm(r0);
iv=v0/norm(v0);
v=norm(v0);
ith=cross(ih,ir);
gamma=atan2(cross(ir,iv),dot(ir,iv));
singamma=norm(cross(ir,iv))/(norm(ir)*norm(iv));
cosgamma=dot(ir,iv);
sinf=(h*v/mu)*(cosgamma/e);
cosf=(1/e)*(h*v*singamma/mu-1);

f0=atan2(sinf,cosf);
if f0<0
    f0=f0+2*pi;
end
if f0>2*pi
    f0=f0-2*pi;
end

sinE=sqrt(1-e^2)*sinf/(1+e*cosf);
cosE=(e+cosf)/(1+e*cosf);
E0=atan2(sinE,cosE);
if E0<0
    E0=E0+2*pi;
end
if E0>2*pi
    E0=E0-2*pi;
end
% E0=acos((e+cos(f0))/(1+e*cos(f0)));
% E0=arc_angles('acos',(e+cos(f0))/(1+e*cos(f0)));
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
