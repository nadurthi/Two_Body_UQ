%% eg2 from doostan paper
r1 = [6777.828 1085.564 -210.326 -0.688773 5.351702 5.387790]';
r2 = [6780.042 1068.401 -227.858 -0.659875 5.357861 5.385113]';
opt = odeset('reltol',1e-12,'abstol',1e-12);
T=[0,172740:1:172860];
% T=0:30:10000;
% [t,X1]= ode45(@twoBody,T,r1,opt);
% [t,X2]= ode45(@twoBody,T,r2,opt);
% [r11,v11]=FGsol(T(end),0,r1(1:3),r1(4:6));
% [r11,X1(end,1:3)']

%% molniya orbit
r1=[-10515.45,-5235.37,49.17,-2.10305,-4.18146,5.563290];
T=0:30:172860;
opt = odeset('reltol',1e-12,'abstol',1e-12);
[t,X1]= ode45(@twoBody,T,r1,opt);
% plot3(X1(:,1),X1(:,2),X1(:,3))
P1=blkdiag(0.01,0.01,0.01,1e-007,1e-007,1e-007);
r2=mvnrnd(r1,P1);
[t,X2]= ode45(@twoBody,T,r2,opt);
plot3(X1(:,1),X1(:,2),X1(:,3),'r',X2(:,1),X2(:,2),X2(:,3),'b')
plot(t,sqrt(sum((X1(:,[1,2,3])-X2(:,[1,2,3])).^2,2)))
%%
% X=zeros(length(T),3);
% i=1;
% for t=T
% [r11,v11]=FGsol(t,0,r1(1:3),r1(4:6));
% 
% X(i,:)=r11';
% i=i+1;
% % plot3(X(:,1),X(:,2),X(:,3))
% % pause(0.1)
% end

% N=1000;
% %T1=0:30:160740;
% %T2=0:30:21600;
% T1=[0:30:172800];
% T2=[0:30:172800];


P1=blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
P2 =blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
% [X01gh5,w1gh5]=GH_points(r1,P1,5);
% [X02gh5,w2gh5]=GH_points(r2,P2,5);
% N=length(w1gh5);
% X1gh5=zeros(length(T),6,N);
% X2gh5=zeros(length(T),6,N);
% parfor i=1:N
%     i
%    [t,X1]= ode45(@twoBody,T,X01gh5(i,:)',opt); 
%    [t,X2]= ode45(@twoBody,T,X02gh5(i,:)',opt);
%    X1gh5(:,:,i)=X1;
%    X2gh5(:,:,i)=X2;
% end
% 
% [X01cut6,w1cut6]=conjugate_dir_gausspts_till_6moment_scheme2(r1,P1);
% [X02cut6,w2cut6]=conjugate_dir_gausspts_till_6moment_scheme2(r2,P2);
% N=length(w1cut6);
% X1cut6=zeros(length(T),6,N);
% X2cut6=zeros(length(T),6,N);
% parfor i=1:N
%     i
%    [t,X1]= ode45(@twoBody,T,X01cut6(i,:)',opt); 
%    [t,X2]= ode45(@twoBody,T,X02cut6(i,:)',opt);
%    X1cut6(:,:,i)=X1;
%    X2cut6(:,:,i)=X2;
% end
% 
% [X01cut8,w1cut8]=conjugate_dir_gausspts_till_8moment(r1,P1);
% [X02cut8,w2cut8]=conjugate_dir_gausspts_till_8moment(r2,P2);
% N=length(w1cut8);
% X1cut8=zeros(length(T),6,N);
% X2cut8=zeros(length(T),6,N);
% parfor i=1:N
%     i
%    [t,X1]= ode45(@twoBody,T,X01cut8(i,:)',opt); 
%    [t,X2]= ode45(@twoBody,T,X02cut8(i,:)',opt);
%    X1cut8(:,:,i)=X1;
%    X2cut8(:,:,i)=X2;
% end
% save('Doostaneg_GH_CUT6_CUT8','X1cut8','w1cut8','X2cut8','w2cut8','X1cut6','w1cut6','X2cut6','w2cut6','X1gh5','w1gh5','X2gh5','w2gh5','T','-v7.3')

% r1 = [6777.828 1085.564 -210.326 -0.688773 5.351702 5.387790]';
% r2 = [6780.042 1068.401 -227.858 -0.659875 5.357861 5.385113]';
% P1
% P2
% X1mc0=mvnrnd(r1,P1,N);
% X2mc0=mvnrnd(r2,P2,N);
% 
% X1mc=zeros(length(T1),6,N);
% X2mc=zeros(length(T2),6,N);
% parfor i=1:N
%     i
%    [t,X1]= ode45(@twoBody,T1,X1mc0(i,:)',opt); 
%    [t,X2]= ode45(@twoBody,T2,X2mc0(i,:)',opt);
%    X1mc(:,:,i)=X1;
%    X2mc(:,:,i)=X2;
% end
% 
   
% N=size(X1mc,3);
% t=size(X1mc,1);
% X1=zeros(N,6);
% X2=zeros(N,6);
% for i=1:1:t
% for j=1:1:N
%         X1(j,:)=X1mc(i,:,j);
%         X2(j,:)=X2mc(i,:,j);
% end 
% %     plot3(X1(:,1),X1(:,2),X1(:,3),'ro',X2(:,1),X2(:,2),X2(:,3),'bo')
% %     grid on
% figure(1)
%     clf
%     subplot(1,3,1)
%     plot(X1(:,1),X1(:,2),'ro',X2(:,1),X2(:,2),'bo')
%     hold on
%     plot_ellipse_cov(mean([X1(:,1),X1(:,2)],1),cov([X1(:,1),X1(:,2)]),[1,2,3,4])    
%     title(num2str(i))
%     
%     subplot(1,3,2)
%     plot(X1(:,2),X1(:,3),'ro',X2(:,2),X2(:,3),'bo')
%     hold on
%     plot_ellipse_cov(mean([X1(:,2),X1(:,3)],1),cov([X1(:,2),X1(:,3)]),[1,2,3,4])
%     title(num2str(i))
%     
%     subplot(1,3,3)
%     plot(X1(:,1),X1(:,3),'ro',X2(:,1),X2(:,3),'bo')
%     hold on
%     plot_ellipse_cov(mean([X1(:,1),X1(:,3)],1),cov([X1(:,1),X1(:,3)]),[1,2,3,4])
%     title(num2str(i))
%     pause(0.1)
% end
% 
% 
