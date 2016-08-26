function [r,v]=FGsol(t,t0,r0,v0)
r0=r0(:);
nr0=norm(r0);
v0=v0(:);
mu=398601.2;
sig0=dot(r0,v0)/sqrt(mu);
a=1/(2/nr0-norm(v0)^2/mu);
h=cross(r0,v0);
e=sqrt(1-norm(h)^2/mu*(2/nr0-norm(v0)^2/mu));
p=norm(h)^2/mu;

sinf0=sqrt(p)*sig0/nr0;
cosf0=p/nr0-1;

f0=atan2(sinf0,cosf0);
if f0<0
    f0=f0+2*pi;
end
opt = odeset('reltol',1e-12,'abstol',1e-12);
c=sqrt(mu/p^3);
if t==t0
    f=f0;
else
[tt,f]=ode45(@(t,f)c*(1+e*cos(f))^2,[t0,t],f0,opt);
f=f(end);
end
if f<0 && f>-2*pi
    f=f+2*pi;
end
if f<-2*pi && f>-4*pi
    f=f+4*pi;
end


th=f-f0;

r=p*nr0/(nr0+(p-nr0)*cos(th)-sqrt(p)*sig0*sin(th));

F=1-r/p*(1-cos(th));
G=r*nr0*sin(th)/sqrt(mu*p);
Ft=sqrt(mu)*(sig0*(1-cos(th))-sqrt(p)*sin(th))/(nr0*p);
Gt=1-nr0/p*(1-cos(th));


r=F*r0+G*v0;
v=Ft*r0+Gt*v0;
end







