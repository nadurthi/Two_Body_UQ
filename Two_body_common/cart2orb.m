function [a,e,i,omg,Omg,M]=cart2orb(x,y,z,xd,yd,zd)

r0=[x;y;z];
v0=[xd;yd;zd];
mu=398600.8;
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

ir=r0/norm(r0);
iv=v0/norm(v0);
v=norm(v0);
singamma=norm(cross(ir,iv))/(norm(ir)*norm(iv));
cosgamma=dot(ir,iv);
sinf=(h*v/mu)*(cosgamma/e);
cosf=(1/e)*(h*v*singamma/mu-1);

sinE=sqrt(1-e^2)*sinf/(1+e*cosf);
cosE=(e+cosf)/(1+e*cosf);
E0=atan2(sinE,cosE);
if E0<0
    E0=E0+2*pi;
end
if E0>2*pi
    E0=E0-2*pi;
end


M=E0-e*sin(E0);


