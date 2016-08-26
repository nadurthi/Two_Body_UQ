% Trying to estimate the pdf in all 6 dimensions in equinoctial element
% space
T=699;
X=zeros(100000,6);
for i=1:1:100000
    X(i,:)=Xmc(T,:,i);
end

X_check=zeros(size(X));
% X_equinoc=zeros(size(X));

            for j=1:1:1000
                S=X_equinoc(j,1:6);
%                   S=X(j,1:6);
%                 [a,P1,P2,Q1,Q2,l]=cart2equinf2(S(1),S(2),S(3),S(4),S(5),S(6));             
                [x,y,z,xd,yd,zd]=equin2cart(S(1),S(2),S(3),S(4),S(5),S(6));
                X_check(j,:)=[x,y,z,xd,yd,zd];
%             X_equinoc(j,:)=[a,P1,P2,Q1,Q2,l];
            end
%  figure
%  plot3(X_equinoc(:,1),X_equinoc(:,2),X_equinoc(:,3),'bo')
%  figure
%  plot3(X_equinoc(:,2),X_equinoc(:,3),X_equinoc(:,4),'bo')
%   figure
%  plot3(X_equinoc(:,3),X_equinoc(:,4),X_equinoc(:,5),'bo')
%  figure
%  plot3(X_equinoc(:,4),X_equinoc(:,5),X_equinoc(:,6),'bo')
%   figure
%  plot3(X_equinoc(:,1),X_equinoc(:,2),X_equinoc(:,4),'bo')
%    figure
%  plot3(X_equinoc(:,1),X_equinoc(:,2),X_equinoc(:,5),'bo')
%    figure
%  plot3(X_equinoc(:,1),X_equinoc(:,2),X_equinoc(:,6),'bo')
%    figure
%  plot3(X_equinoc(:,2),X_equinoc(:,3),X_equinoc(:,5),'bo')
%     figure
%  plot3(X_equinoc(:,2),X_equinoc(:,3),X_equinoc(:,6),'bo')
%     figure
%  plot3(X_equinoc(:,3),X_equinoc(:,4),X_equinoc(:,6),'bo')
%     figure
%  plot3(X_equinoc(:,2),X_equinoc(:,3),X_equinoc(:,5),'bo'
% Xxy=X(1,:); 
clc
P=0;
for j=1:1:10000
%     j
Xxy=Xmc_1e4(1,:,j) ;
% Xxy=[7000 0 0 0.1 -1.0374090357 7.4771288355];
%  Xxy=[1000,0,0,0,sqrt(398601.2/1000)+4,0];
 [a,P1,P2,Q1,Q2,l]=cart2equin(Xxy(1),Xxy(2),Xxy(3),Xxy(4),Xxy(5),Xxy(6));
[x,y,z,xd,yd,zd]=equin2cart(a,P1,P2,Q1,Q2,l);
S=max(abs( Xxy-[x,y,z,xd,yd,zd]));
if S>=P
    S
    P=S;
end
end

%  Xxy=[1000,0,0,0,sqrt(398601.2/1000)+4,0];
   opt = odeset('reltol',1e-10,'abstol',1e-10);
 [t,Xplot]= ode45(@twoBody,[0 20000],Xxy',opt);
 plot3(Xplot(:,1),Xplot(:,2),Xplot(:,3))
%  axis square
%  
[x,y,z,xd,yd,zd]=equin2cart(a,P1,P2,Q1,Q2,l);
 [an,P1n,P2n,Q1n,Q2n,ln]=cart2equin(x,y,z,xd,yd,zd);
max(abs([an,P1n,P2n,Q1n,Q2n,ln]-[a,P1,P2,Q1,Q2,l]))
