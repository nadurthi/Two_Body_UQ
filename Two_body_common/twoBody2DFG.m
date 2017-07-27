
function x = twoBody2DFG(t,t0,x)
r0=x(1:2);
v0=x(3:4);
r0=r0(:)';
v0=v0(:)';

options=optimset('TolCon',1e-12,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-15,'TolX',1e-15);



mu=398601.2; %km^3/sec^2
a=(2/norm(r0)-norm(v0)^2/mu)^(-1);
h=cross([r0,0],[v0,0]);
ve=1/mu*(cross([v0,0],h)-mu/norm(r0)*[r0,0]);
e=norm(ve);
sig0=dot([r0,0],[v0,0])/sqrt(mu);


Re=6378.165; %km

Ec=0;
Ec=fsolve(@(E)[sqrt(mu/a^3)*(t-t0)-E+e*sin(E)],Ec,options);

%%%%%%%%%%%%%%% F anf G solution
    r=a+(norm(r0)-a)*cos(Ec)+sqrt(a)*sig0*sin(Ec);
%     r_norm(i)=r;
    F=1-a/norm(r0)*(1-cos(Ec));
    Ft=-sqrt(mu*a)/(r*norm(r0))*sin(Ec);
%     G=a*sig0/sqrt(mu)*(1-cos(Ec))+norm(r0v)*sqrt(a/mu)*sin(Ec);
    G=(t-t0)+sqrt(a^3/mu)*(sin(Ec)-Ec);
    Gt=1-a/r*(1-cos(Ec));
    x(1:2)=(F*r0+G*v0)';
    x(3:4)=(Ft*r0+Gt*v0)';
    
end