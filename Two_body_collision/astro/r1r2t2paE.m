%  [p,a,dE,v11,v22] = r1r2t2paE(r1v,r2v,tx,mu,shortorlong)
%
%   uses the p-iteration method
%
%			r1v, r2v are the two vectors
%			tx is the time of flight
%			mu is the two-body gravitational constant
%			shortorlong = 0|1 (short|long)
%		

function [p,a,dE,v11,v22] = r1r2t2paE(r1v,r2v,tx,mu,shortorlong)
	r1   = norm(r1v);
	r2   = norm(r2v);
   dnu  = acos( (r1v'*r2v)/r1/r2 );
   if nargin == 5,
		if shortorlong == 1,
			dnu = 2*pi-dnu;
		end;
	end;
	k    = r1*r2*(1-cos(dnu));
	l    = r1+r2;
	m    = r1*r2*(1+cos(dnu));
	ppi  = k/(l+sqrt(2*m));
	ppii = k/(l-sqrt(2*m));
	po   = (2*ppi +   ppii)/3;
   p1   = (  ppi + 2*ppii)/3;
   
	ao   = m*k*po/( (2*m-l*l)*po^2 + 2*k*l*po - k^2);
	fo   = 1 - r2/po *(1-cos(dnu));
   dEo  = atan2(-r1*r2/sqrt(mu*ao),1-r1/ao*(1-fo));
   if (dEo < 0), dEo = dEo+2*pi; end
	go   = r1*r2*sin(dnu)/sqrt(mu*po);
	to   = go + sqrt(ao^3/mu)*(dEo-sin(dEo));

	a1   = m*k*p1/( (2*m-l*l)*p1^2 + 2*k*l*p1 - k^2);
	f1   = 1 - r2/p1 *(1-cos(dnu));
   dE1  = atan2(-r1*r2/sqrt(mu*a1),1-r1/a1*(1-f1));
   if (dE1 < 0), dE1 = dE1+2*pi; end
	g1   = r1*r2*sin(dnu)/sqrt(mu*p1);
   t1   = g1 + sqrt(a1^3/mu)*(dE1-sin(dE1));
   
   %figure(1); hdl=plot(po,ao,'o'); set(hdl,'markersize',20);
   %figure(1); hdl=plot(p1,a1,'o'); set(hdl,'markersize',20);

   %figure(2); hdl=plot(p1,t1,'o'); set(hdl,'markersize',20);

   dt   = t1 - tx;
   
   its = 0;  maxits = 100;
   savem = [p1 a1 dt];
   while ( abs(dt) > 1e-11 & its < maxits )
      its = its+1;
      deltap = - dt*(p1-po)/(t1-to);
		p = p1 + deltap/(1+deltap^2);
		po = p1; to = t1;
		p1 = p;
		a1   = m*k*p1/( (2*m-l*l)*p1^2 + 2*k*l*p1 - k^2);
      f1   = 1 - r2/p1 *(1-cos(dnu));
      fd1  = sqrt(mu/p1)*tan(dnu/2)*((1-cos(dnu))/p1-1/r1-1/r2);
		g1   = r1*r2*sin(dnu)/sqrt(mu*p1);
      if ( a1 > 0 )
         dE1  = atan2(-r1*r2/sqrt(mu*a1)*fd1,1-r1/a1*(1-f1));
         if (dE1<0), dE1=dE1+2*pi;end
         t1   = g1 + sqrt(a1^3/mu)*(dE1-sin(dE1));
      else
         dF1 = acosh(1-r1/a1*(1-f1));
         t1  = g1 + sqrt((-a1)^3/mu)*(sinh(dF1)-dF1);
      end
      dt = t1-tx;
      dtshow   = [t1  tx t1-tx];
      %figure(1);plot(p,a1,'o');
      %figure(2); plot(p,t1,'o');
      %drawnow;pause
      savem = [savem; p a1 dt];
   end;
   
   if (a1 > 0 )
      a = a1; dE = dE1;
   else
      a = a1; dE = dF1;
   end
   %savem
   f1;
   fd1;
   g1;
   gdot=(1+fd1*g1)/f1;
   v11 = (r2v-f1*r1v)/g1;
   v22=(fd1-f1*gdot/g1)*r1v+gdot/g1*r2v;
   %r22=f1*r1v+g1*v11
   %r1*r2*(1-cos(dnu))/a/p/(1-cos(dE))