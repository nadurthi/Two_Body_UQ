function [Omg,omg,i]=solve_rot_mat(ix,iy,iz,ie,ip,ih)
l1=dot(ix,ie);
l2=dot(ix,ip);
l3=dot(ix,ih);

m1=dot(iy,ie);
m2=dot(iy,ip);
m3=dot(iy,ih);

n1=dot(iz,ie);
n2=dot(iz,ip);
n3=dot(iz,ih);
R=[l1,l2,l3;m1,m2,m3;n1,n2,n3];
% 
% omg=atan2(n1,n2);
% 
% si1=sqrt(n1^2+n2^2);
% si2=-sqrt(n1^2+n2^2);
% 
% ci=n3;
% i1=atan2(si1,ci);
% i2=atan2(si2,ci);
% 
% Omg1=atan2(-l3,m3);
% Omg2=atan2(l3,-m3);
% 
% si3=sqrt(l3^2+m3^2);
% si4=-sqrt(l3^2+m3^2);
% 
% i3=atan2(si3,ci);
% i4=atan2(si4,ci);
% 
% i=[i1,i2,i3,i4];
% Omg=[Omg1,Omg2];
% f=zeros(4,2);
% mini=10000;
% minOmg=10000;
% f1=100000000;
% 
% for ii=1:1:4
%     for j=1:2
%     e1=abs(l1-(cos(Omg(j))*cos(omg)-sin(Omg(j))*sin(omg)*cos(i(ii)))); 
%     e2=abs(l2-(-cos(Omg(j))*sin(omg)-sin(Omg(j))*cos(omg)*cos(i(ii))));
%     e3=abs(l3-(sin(Omg(j))*sin(i(ii))));
%     e4=abs(m1-(sin(Omg(j))*cos(omg)+cos(Omg(j))*sin(omg)*cos(i(ii))));
%     e5=abs(m2-(-sin(Omg(j))*sin(omg)+cos(Omg(j))*cos(omg)*cos(i(ii))));
%     e6=abs(m3-(-cos(Omg(j))*sin(i(ii))));
%     e7=abs(n1-(sin(omg)*sin(i(ii))));
%     e8=abs(n2-(cos(omg)*sin(i(ii))));
%     e9=abs(n3-cos(i(ii)));
%     f(ii,j)=e1+e2+e3+e4+e5+e6+e7+e8+e9;
%     if f(ii,j)<=f1
%         f1=f(ii,j);
%     mini=ii;
%     minOmg=j;
%     end
%     
%     end
% end
% i=i(mini);
% Omg=Omg(minOmg);
% 

OUTPUT=SpinCalc('DCMtoEA313',R',1e-13,0)*pi/180;
Omg=OUTPUT(1);
omg=OUTPUT(3);
i=OUTPUT(2);
if i<0
    i=i+2*pi;
end
if i>2*pi
    i=i-2*pi;
end    
if Omg<0
    Omg=Omg+2*pi;
end
if Omg>2*pi
    Omg=Omg-2*pi;
end
if omg<0
    omg=omg+2*pi;
end
if omg>2*pi
    omg=omg-2*pi;
end
% max(abs([Omg-Omg2,i-i2,omg-omg2]))
end