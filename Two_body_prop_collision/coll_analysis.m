% creating a new problem for collision
P0 =blkdiag(0.01,0.01,0.01,0.0000001,0.0000001,0.0000001);
% r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
% r2 = [7000 0 0 0 7.4771288355  1.0374090357]';
r1=[6777.828,1085.564,-210.326,-0.688773,5.351702,5.387790];
r2=[6780.042,1068.401,-227.858,-0.659875,5.357861,5.385113];

opt = odeset('reltol',1e-15,'abstol',1e-15);

time.t0 = 0;
time.dt = 30;
time.tf = 172800+300;
% T=[linspace(time.t0,3*time.tf/4,200),linspace(3*time.tf/4+time.dt,time.tf,400)];
T=time.t0:time.dt:time.tf;
Tb=time.tf:-time.dt:time.t0;
% 
% [t,X1]= ode45(@twoBody,Tb,r1,opt);
% [t,X2]= ode45(@twoBody,Tb,r2,opt);
% 
% d=zeros(length(T),1);
% for i=1:1:length(Tb)
% 
%     plot3(X1(1:i,1),X1(1:i,2),X1(1:i,3),'ro',X2(1:i,1),X2(1:i,2),X2(1:i,3),'bo')
%     pause(0.01)
%     d(i)=norm([X1(i,1),X1(i,2),X1(i,3)]-[X2(i,1),X2(i,2),X2(i,3)]);
% end
% 
% rr1=[X1(end,1:3),X1(end,4:6)];
% rr2=[X2(end,1:3),X2(end,4:6)];

[t,X11]= ode45(@twoBody,T,r1,opt);
[t,X22]= ode45(@twoBody,T,r2,opt);

d=zeros(length(T),1);
for i=1:1:length(T)

%     plot3(X11(1:i,1),X11(1:i,2),X11(1:i,3),'ro',X22(1:i,1),X22(1:i,2),X22(1:i,3),'bo')
%     pause(0.01)
    d(i)=norm([X11(i,1),X11(i,2),X11(i,3)]-[X22(i,1),X22(i,2),X22(i,3)]);
end
rrtr1=rr1'+sqrtm(P0)*randn(6,1);
rrtr2=rr2'+sqrtm(P0)*randn(6,1);

%%
N=5000;
wmc=ones(N,1)./N;
X1mc0=mvnrnd(rrtr1,P0,N);
X2mc0=mvnrnd(rrtr2,P0,N);

X1mc=zeros(length(T),6,N);
X2mc=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,X1mc0(i,:)',opt);
   X1mc(:,:,i)=X;
   [t,X]= ode45(@twoBody,T,X2mc0(i,:)',opt);
   X2mc(:,:,i)=X;
end
save('Collision2sims','X1mc','X2mc','wmc')

X1=zeros(N,6);
X2=zeros(N,6);
for i=570:1:599
    for j=1:1:N
        X1(j,:)=X1mc(i,:,j);
        X2(j,:)=X2mc(i,:,j);
    end
   plot3(X1(:,1),X1(:,2),X1(:,3),'ro',X2(:,1),X2(:,2),X2(:,3),'bo') 
%    axis([-1e4,1e4,-1e4,1e4,-1e4,1e4])
    pause(1)
end

d=zeros(N,1);
D=zeros(N,3);
ind=randperm(N);
% k=1;
n=0;
for i=1:1:N
%     for j=1:1:N
    d(i)=norm(X1(ind(i),1:3)-X2(i,1:3));
    D(i,:)=(X1(ind(i),1:3)-X2(i,1:3))/1;
    if d(i)<1
        n=n+1;
    end
%     k=k+1;
%     end
end
 hist(d,100)
pd = fitdist(d,'Gamma');
x=0:0.01:5;
pf = pdf(pd,x);
figure
plot(x,pf)
figure
plot3(D(:,1),D(:,2),D(:,3),'ro')

