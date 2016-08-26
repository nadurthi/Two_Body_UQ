function ZZ = XYZ2OE_multiple(YY)
N=size(YY,1);
ZZ=zeros(size(YY));
for kl=1:1:N
    XX=YY(kl,:);
XX=XX(:);
X=XX(1:3);
V=XX(4:6);

 mu=398601.2;
 
%% Input - X,V position and velocity of the point mass... 

rm = sqrt(X'*X); vm = sqrt(V'*V);

a = 1/(2/rm - vm^2/mu);

hbar = cross(X,V);
cbar = cross(V,hbar)-mu*X/rm;

e = sqrt(cbar'*cbar)/mu;

hm = sqrt(hbar'*hbar); ih = hbar/hm;

ie = cbar/(mu*e);

ip = cross(ih,ie);

i = acos(ih(3,1));

w = atan2(ie(3,1),ip(3,1));

Om = atan2(ih(1,1),-ih(2,1));

sig = X'*V/sqrt(mu);

E = atan2(sig/sqrt(a),1-rm/a);

OE=[a,e,E,w,i,Om];
ZZ(kl,:)=OE;
end
