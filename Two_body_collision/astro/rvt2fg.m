%  rvt2fg.m  computes f,g,fdot,gdot from initial r,v,dt
%
%	[fg,r,v] = rvt2fg(ro,vo,dt,mu)
%		ro,vo = initial column vectors in IJK frame
%		dt = t-to
%		fg = [f g fdot gdot]
%
function [fg,r,v] = rvt2fg(ro,vo,dt,mu)
	ros= sqrt(ro'*ro);
	aen = rv2aen(ro,vo,mu);
	a  = aen(1);  e = aen(2);  nuo = aen(3);
	p  = a*(1-e^2);
	n  = sqrt(mu/a^3);
    PP = 2*pi*sqrt(a^3/mu);
    while ( dt > PP )
        dt = dt - PP;
    end
    dM = n*dt;
   
   Eo = 2*atan(sqrt((1-e)/(1+e))*tan(nuo/2)) ;
   M=dM+Eo-e*sin(Eo);
   E1 = EofMe(M,e,1.0e-11);
   rs = a*(1-e*cos(E1));
	dE = E1-Eo ;
	f  = 1 - a/ros *(1-cos(dE));
	g  = dt - (dE-sin(dE))/n;
	fd = -sqrt(mu*a)*sin(dE)/rs/ros;
	gd = 1 - a/rs *(1-cos(dE));
	fg = [f g fd gd];
	r  = f*ro +  g*vo;
	v  = fd*ro + gd*vo;
