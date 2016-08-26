
% 2/09/2009 at 11:57:36.8904960001419
CML1='1 22729U 93049A   13365.23601771  .99999999  96287-5  94692-3 0  8555';
CML2='2 22729 062.1029 006.1492 0240560 285.6452 156.4968 16.15760818150519     0.00      4320.0        360.00';
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
   global opsmode
   opsmode= 'a';
global idebug dbgfile

endtime.year=2014;
endtime.month=1;
endtime.day=15;
endtime.hr=1;
endtime.min=1;
endtime.sec=1;
dtmin=0.5;
[satrecCM, startmfeCM, stopmfeCM, dtmin] = twoline2rv_modified(72,CML1,CML2,'m','e',endtime,dtmin);
 
 [satrecCM, rCM0, vCM0] = sgp4(satrecCM,0);


 T=0:0.1:500;
XtestCM=zeros(length(T),6);
 for t =1:1:length(T)
    
    [satrecCMmc, rCM, vCM] = sgp4(satrecCM,stopmfeCM+T(t));
    XtestCM(t,:)=[rCM(:)',vCM(:)'];
    
 end
plot3(XtestCM(:,1),XtestCM(:,2),XtestCM(:,3),'ro-')
hold on
[x,y,z] = sphere;
surf(x*6378.137,y*6378.137,z*6378.137) 


% r0=[-10515.45,-5235.37,49.17];
% v0=[-2.10305,-4.18146,5.563290];
r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
r2 = [7000 0 0 0 7.4771288355  1.0374090357]';
opt = odeset('reltol',1e-12,'abstol',1e-12);
[t,X11]= ode45(@twoBody,T*60,[r0,v0],opt);
for i=1:10:length(T)
plot3(X11(1:i,1),X11(1:i,2),X11(1:i,3),'ro-')
title(num2str(i))
hold on
surf(x*6378.137,y*6378.137,z*6378.137) 
hold off
pause(0.01);
end

N=5000;
 T=0:30:172800;
P0 =blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
X1mc0=mvnrnd(r1,P0,N);
X1mc=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,X1mc0(i,:)',opt);
   X1mc(:,:,i)=X;
end

N=size(X1mc,3);
X1=zeros(N,6);
tt=[2695:2730,2895:2930,3095:3030,3295:3330,3395:3430,3595:3630,3795:3830,3995:4030,4195:4230,5350:5370,5500:5540,5700:5730];
for i=tt
    for j=1:1:N
        X1(j,:)=X1mc(i,:,j);
    end
    figure(1)
    clf
    subplot(1,3,1)
    plot(X1(:,1),X1(:,2),'ro')
    hold on
    plot_ellipse_cov(mean([X1(:,1),X1(:,2)],1),cov([X1(:,1),X1(:,2)]),[1,2,3,4])    
    title(num2str(i))
    
    subplot(1,3,2)
    plot(X1(:,2),X1(:,3),'ro')
    hold on
    plot_ellipse_cov(mean([X1(:,2),X1(:,3)],1),cov([X1(:,2),X1(:,3)]),[1,2,3,4])
    title(num2str(i))
    
    subplot(1,3,3)
    plot(X1(:,1),X1(:,3),'ro')
    hold on
    plot_ellipse_cov(mean([X1(:,1),X1(:,3)],1),cov([X1(:,1),X1(:,3)]),[1,2,3,4])
    title(num2str(i))
    saveas(gca,strcat('firstbody_',num2str(i)),'jpg')
pause(0.1)
end

%% second body
rr2=[197,853.4,-6097,8.56477302829007,0.259513978551832, 0.0561224124116941];
% rr2=[35.8,849,-6097,0.0332421258473425,-1.17930717307689, 8.4876782282537];
 T=0:30:(5359*30);
 TT=(5359*30):-30:0;
 [t,X]= ode45(@twoBody,TT,rr2,opt);
 r2=X(end,:);
 [t,X]= ode45(@twoBody,T,r2,opt);
 
 
 N=100;
P0 =blkdiag(0.01,0.01,0.01,1e-0010,1e-0010,1e-0010);
X2mc0=mvnrnd(r2,P0,N);
X2mc0=[r2(:)';X2mc0];
X2mc0(end,:)=[];

X2mc=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,X2mc0(i,:)',opt);
   X2mc(:,:,i)=X;
end
 
for j=1:1:N
        X2(j,:)=X2mc(end,:,j);
end

plot(X1(:,1),X1(:,3),'ro')
    hold on
plot_ellipse_cov(mean([X1(:,1),X1(:,3)],1),cov([X1(:,1),X1(:,3)]),[1,2,3,4])
title(num2str(i))
plot(X2(:,1),X2(:,3),'bo')
plot_ellipse_cov(mean([X2(:,1),X2(:,3)],1),cov([X2(:,1),X2(:,3)]),[1,2,3,4])

d=zeros(N,1);
D=zeros(N,3);
ind=randperm(N);
% k=1;
n=0;
for i=1:1:N
%     for j=1:1:N
    d(i)=norm(X1(ind(i),1:3)-X2(i,1:3));
    D(i,:)=(X1(ind(i),1:3)-X2(i,1:3))/1;
    if d(i)<0.5
        n=n+1;
    end
%     k=k+1;
%     end
end

% gauss coll
% XG1=mvnrnd(mean(X1(:,1:3)),cov(X1(:,1:3)),5000);
% XG2=mvnrnd(mean(X2(:,1:3)),cov(X2(:,1:3)),5000);

dg=zeros(N,1);
Dg=zeros(N,3);
ind=randperm(N);
% k=1;
ng=0;
for i=1:1:N
%     for j=1:1:N
    dg(i)=norm(XG1(ind(i),1:3)-XG2(i,1:3));
    Dg(i,:)=(XG1(ind(i),1:3)-XG2(i,1:3))/1;
    if dg(i)<0.5
        ng=ng+1;
    end
%     k=k+1;
%     end
end
[n,ng]
figure
hist(dg)
figure
hist(d)