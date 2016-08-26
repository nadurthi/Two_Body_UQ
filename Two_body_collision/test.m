% r0=[4743.71416928;3008.08450440;843.2102315222];
% v0=[-4.7153;6.0107;5.0842];
% [a,P1,P2,Q1,Q2,l]=cart2equin(4743.71416928,3008.08450440,843.2102315222,-4.7153,6.0107,5.0842)
% [x,y,z,xd,yd,zd]=equin2cart(a,P1,P2,Q1,Q2,l)
% [a,e,i,omg,Omg,M]=equin2orbital(a,P1,P2,Q1,Q2,l)
% [a,P1,P2,Q1,Q2,l]=orbital2equin(a,e,i,omg,Omg,M)


r0 = [6999.97317756185,0.187428082420522,0.12817388014082,0.001010925130436,-1.03621859775461,7.47932400986844]%[7000 0 0 0 -1.0374090357 7.4771288355]';
[a,P1,P2,Q1,Q2,l]=cart2equin(r0(1),r0(2),r0(3),r0(4),r0(5),r0(6));

[x,y,z,xd,yd,zd]=equin2cart(a,P1,P2,Q1,Q2,l);
 [an,P1n,P2n,Q1n,Q2n,ln]=cart2equin(x,y,z,xd,yd,zd);
max(abs([an,P1n,P2n,Q1n,Q2n,ln]-[a,P1,P2,Q1,Q2,l]))

time.t0 = 0;
time.dt = 29.173785213265145;
time.tf = 21600;
T=time.t0:time.dt:time.tf;
 opt = odeset('reltol',1e-10,'abstol',1e-10);
 [t,Xxy]= ode45(@twoBody,T,r0,opt);
[t,Xeq]= ode45(@diff_equinoctial,T,[a,P1,P2,Q1,Q2,l]',opt);

Xxytrans=zeros(size(Xeq));
for i=1:1:size(Xeq,1)
    [x,y,z,xd,yd,zd]=equin2cart(Xeq(i,1),Xeq(i,2),Xeq(i,3),Xeq(i,4),Xeq(i,5),Xeq(i,6));
%     [x,y,z,xd,yd,zd]=cart2equin(Xxy(i,1),Xxy(i,2),Xxy(i,3),Xxy(i,4),Xxy(i,5),Xxy(i,6));
    Xxytrans(i,:)=[x,y,z,xd,yd,zd];
end
for i=1:1:length(t)
plot3(Xxy(i:i+10,1),Xxy(i:i+10,2),Xxy(i:i+10,3),'rs-',Xxytrans(i:i+10,1),Xxytrans(i:i+10,2),Xxytrans(i:i+10,3),'bo-','linewidth',2)
% plot3(Xeq(1:i,1),Xeq(1:i,2),Xeq(1:i,3),'r',Xxytrans(1:i,1),Xxytrans(1:i,2),Xxytrans(1:i,3),'bo-','linewidth',2)
axis([ -10000       10000       -1000        1000      -10000       10000])
pause(0.2)
end
figure
plot3(Xxy(:,1),Xxy(:,2),Xxy(:,3),'rs-','linewidth',2,'Markersize',6)
hold on
plot3(Xxytrans(:,1),Xxytrans(:,2),Xxytrans(:,3),'bo-','linewidth',2)
grid
 max(max(Xxy-Xxytrans))

