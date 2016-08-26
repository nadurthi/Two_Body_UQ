% r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
% opt = odeset('reltol',1e-12,'abstol',1e-12);
% T=0:30:172800;
% [t,X1traj]= ode45(@twoBody,T,r1,opt);
% rt=X1traj(5359,1:3);
% 
% vt=X1traj(5359,4:6);
% h=cross(rt,vt);
% h=h/norm(h);
% vt2=0.1*norm(vt)*h+0.98*vt;
%  rt=[189,851.7,-6097.07249833529];
% r2=[rt,vt2];
% T2=0:30:21600;
% Tb=T2(end:-1:1);
% [t,X2rev]= ode45(@twoBody,Tb,r2,opt);
% r20=X2rev(end,:);
% [t,X2for]= ode45(@twoBody,T2,r20,opt);
% plot3(X2for(:,1),X2for(:,2),X2for(:,3),'ro-',X1traj(1:5359,1),X1traj(1:5359,2),X1traj(1:5359,3),'bo-')
% 
% load('F:\2_body_problem\CollisionAnalysis\collision test case 1\secondbody')
N2=size(X2mc,3);
N=size(X1mc,3);
X1=zeros(N,6);
X2=zeros(N2,6);
for j=1:1:N
        X1(j,:)=X1mc(end,:,j);
end    
for j=1:1:N2
        X2(j,:)=X2mc(end,:,j);
end
    figure(1)
    clf
    subplot(1,3,1)
    plot(X1(:,1),X1(:,2),'ro',X2(:,1),X2(:,2),'bo')
    hold on
    plot_ellipse_cov(mean([X1(:,1),X1(:,2)],1),cov([X1(:,1),X1(:,2)]),[1,2,3,4])    
%     title(num2str(i))
    
    subplot(1,3,2)
    plot(X1(:,2),X1(:,3),'ro',X2(:,2),X2(:,3),'bo')
    hold on
    plot_ellipse_cov(mean([X1(:,2),X1(:,3)],1),cov([X1(:,2),X1(:,3)]),[1,2,3,4])
%     title(num2str(i))
    
    subplot(1,3,3)
    plot(X1(:,1),X1(:,3),'ro',X2(:,1),X2(:,3),'bo')
    hold on
    plot_ellipse_cov(mean([X1(:,1),X1(:,3)],1),cov([X1(:,1),X1(:,3)]),[1,2,3,4])
%     title(num2str(i))



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

% gauss coll
XG1=mvnrnd(mean(X1(:,1:3)),cov(X1(:,1:3)),5000);
XG2=mvnrnd(mean(X2(:,1:3)),cov(X2(:,1:3)),5000);

dg=zeros(N,1);
Dg=zeros(N,3);
% ind=randperm(N);
% k=1;
ng=0;
for i=1:1:N
%     for j=1:1:N
    dg(i)=norm(XG1(ind(i),1:3)-XG2(i,1:3));
    Dg(i,:)=(XG1(ind(i),1:3)-XG2(i,1:3))/1;
    if dg(i)<1
        ng=ng+1;
    end
%     k=k+1;
%     end
end
[n,ng]


figure(2)

    clf
    subplot(1,3,1)
    plot(X1(:,1),X1(:,2),'ro',X2(:,1),X2(:,2),'bo',XG1(:,1),XG1(:,2),'+',XG2(:,1),XG2(:,2),'*')
    hold on
    plot_ellipse_cov(mean([X1(:,1),X1(:,2)],1),cov([X1(:,1),X1(:,2)]),[1,2,3,4])    
%     title(num2str(i))
    
    subplot(1,3,2)
    plot(X1(:,2),X1(:,3),'ro',X2(:,2),X2(:,3),'bo',XG1(:,2),XG1(:,3),'+',XG2(:,2),XG2(:,3),'*')
    hold on
    plot_ellipse_cov(mean([X1(:,2),X1(:,3)],1),cov([X1(:,2),X1(:,3)]),[1,2,3,4])
%     title(num2str(i))
    
    subplot(1,3,3)
    plot(X1(:,1),X1(:,3),'ro',X2(:,1),X2(:,3),'bo',XG1(:,1),XG1(:,3),'+',XG2(:,1),XG2(:,3),'*')
    hold on
    plot_ellipse_cov(mean([X1(:,1),X1(:,3)],1),cov([X1(:,1),X1(:,3)]),[1,2,3,4])
%     title(num2str(i))
