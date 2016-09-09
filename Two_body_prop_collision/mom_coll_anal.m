%% moment approach to collision


%% MC aaproach
load('F:\2_body_problem\CollisionAnalysis\CollisionBodies')
N2=size(X2mc,3);
N1=size(X1mc,3);
X1_mc=zeros(N1,3);
X2_mc=zeros(N2,3);
for j=1:1:N1
        X1_mc(j,:)=X1mc(end,1:3,j)/6000;
end    
for j=1:1:N2
        X2_mc(j,:)=X2mc(end,1:3,j)/6000;
end
    figure(1)
    clf
    subplot(1,3,1)
    plot(X1_mc(:,1),X1_mc(:,2),'ro',X2_mc(:,1),X2_mc(:,2),'bo',X1_gh(:,1),X1_gh(:,2),'k+')
    hold on
    plot_ellipse_cov(mean([X1_mc(:,1),X1_mc(:,2)],1),cov([X1_mc(:,1),X1_mc(:,2)]),[1,2,3,4])    
%     title(num2str(i))
    
    subplot(1,3,2)
    plot(X1_mc(:,2),X1_mc(:,3),'ro',X2_mc(:,2),X2_mc(:,3),'bo',X1_gh(:,2),X1_gh(:,3),'k+')
    hold on
    plot_ellipse_cov(mean([X1_mc(:,2),X1_mc(:,3)],1),cov([X1_mc(:,2),X1_mc(:,3)]),[1,2,3,4])
%     title(num2str(i))
    
    subplot(1,3,3)
    plot(X1_mc(:,1),X1_mc(:,3),'ro',X2_mc(:,1),X2_mc(:,3),'bo',X1_gh(:,1),X1_gh(:,3),'k+')
    hold on
    plot_ellipse_cov(mean([X1_mc(:,1),X1_mc(:,3)],1),cov([X1_mc(:,1),X1_mc(:,3)]),[1,2,3,4])
%     title(num2str(i))



d=zeros(N1,1);
D=zeros(N1,3);
ind=randperm(N1);
% k=1;
n=0;
for i=1:1:N1
%     for j=1:1:N
    d(i)=norm(X1_mc(ind(i),1:3)-X2_mc(i,1:3));
    D(i,:)=(X1_mc(ind(i),1:3)-X2_mc(i,1:3))/1;
    if d(i)<0.5/6000
        n=n+1;
    end
%     k=k+1;
%     end
end
XG1=mvnrnd(mean(X1_mc(:,1:3)),cov(X1_mc(:,1:3)),N1);
XG2=mvnrnd(mean(X2_mc(:,1:3)),cov(X2_mc(:,1:3)),N2);

dg=zeros(N1,1);
Dg=zeros(N1,3);
% ind=randperm(N);
% k=1;
ng=0;
for i=1:1:N1
%     for j=1:1:N
    dg(i)=norm(XG1(ind(i),1:3)-XG2(i,1:3));
    Dg(i,:)=(XG1(ind(i),1:3)-XG2(i,1:3))/1;
    if dg(i)<0.5/6000
        ng=ng+1;
    end
%     k=k+1;
%     end
end
[n,ng]

w=ones(N1,1)/N1;

 
 [y1,M1D]=Cal_moments_samples(D,w,1,'raw');
 [y2,M2D]=Cal_moments_samples(D,w,2,'raw');
 [y3,M3D]=Cal_moments_samples(D,w,3,'raw');
 [y4,M4D]=Cal_moments_samples(D,w,4,'raw');

 [y1,M11]=Cal_moments_samples(X1_mc(:,[1:3]),w,1,'raw');
 [y2,M21]=Cal_moments_samples(X1_mc(:,[1:3]),w,2,'raw');
 [y3,M31]=Cal_moments_samples(X1_mc(:,[1:3]),w,3,'raw');
 [y4,M41]=Cal_moments_samples(X1_mc(:,[1:3]),w,4,'raw');

 [y1,M12]=Cal_moments_samples(X2_mc(:,[1:3]),w,1,'raw');
 [y2,M22]=Cal_moments_samples(X2_mc(:,[1:3]),w,2,'raw');
 [y3,M32]=Cal_moments_samples(X2_mc(:,[1:3]),w,3,'raw');
 [y4,M42]=Cal_moments_samples(X2_mc(:,[1:3]),w,4,'raw');
 
%% Gasssian Quadrature Apporach
r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
r2=[6729.43097302094	-318.159560024553	-1381.89229666774	2.26490482713380	-1.18931778372278	7.16940625613625]';                  
% opt = odeset('reltol',1e-12,'abstol',1e-12);
% 
% T1=linspace(0,160740,10);
% T2=linspace(0,21600,10);
% 
% P1=blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
% P2 =blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
% [X1gh0,w1gh] = GH_points(r1,P1,4);
% [X2gh0,w2gh] = GH_points(r2,P2,4);
% 
% X1gh=zeros(length(T1),6,length(w1gh));
% X2gh=zeros(length(T2),6,length(w2gh));
% 
% parfor i=1:length(w1gh)
%     i
%    [t,X]= ode45(@twoBody,T1,X1gh0(i,:)',opt); 
%    X1gh(:,:,i)=X;
%    [t,X]= ode45(@twoBody,T2,X2gh0(i,:)',opt); 
%    X2gh(:,:,i)=X;
% end
% save('Collanal_GH4','X1gh','w1gh','X2gh','w2gh')
load('Collanal_GH4')
X1_gh=zeros(length(w1gh),3);
X2_gh=zeros(length(w1gh),3);
for j=1:1:length(w1gh)
        X1_gh(j,:)=X1gh(end,1:3,j)/6000;
        X2_gh(j,:)=X2gh(end,1:3,j)/6000;
end

 [y1,M11_gh]=Cal_moments_samples(X1_gh(:,[1:3]),w1gh,1,'raw');
 [y2,M21_gh]=Cal_moments_samples(X1_gh(:,[1:3]),w1gh,2,'raw');
 [y3,M31_gh]=Cal_moments_samples(X1_gh(:,[1:3]),w1gh,3,'raw');
 [y4,M41_gh]=Cal_moments_samples(X1_gh(:,[1:3]),w1gh,4,'raw');

 [y1,M12_gh]=Cal_moments_samples(X2_gh(:,[1:3]),w2gh,1,'raw');
 [y2,M22_gh]=Cal_moments_samples(X2_gh(:,[1:3]),w2gh,2,'raw');
 [y3,M32_gh]=Cal_moments_samples(X2_gh(:,[1:3]),w2gh,3,'raw');
 [y4,M42_gh]=Cal_moments_samples(X2_gh(:,[1:3]),w2gh,4,'raw');