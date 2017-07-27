%2 Body problem entropy reconstruction

% clear
% close all
% clc

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('body2sims')
DX=Xut(:,[2,3],:);
w=wut;

Tp=699;
pts=Xmc;
%%%%%%%%%%%%%%%%%%%%%%%%%%

ns=size(DX,2);
N=size(DX,3);

XY=zeros(N,ns);
%  XYpts=zeros(1e5,ns);
 for i=1:1:N
 XY(i,1)=DX(Tp,1,i);
 XY(i,2)=DX(Tp,2,i);
%  XY(i,3)=DX(Tp,3,i);
 end 
% 
%   for i=1:1:100000
%   XYpts(i,1)=Xmc(Tp,1,i);
%  XYpts(i,2)=Xmc(Tp,2,i);
%  XYpts(i,3)=Xmc(Tp,3,i);
%  
%  end
%    XYmc=transform_domain(XYpts(:,1:2),bl,bu,-s*ones(1,ns),s*ones(1,ns));
%    YZmc=transform_domain(XYpts(:,2:3),bl,bu,-s*ones(1,ns),s*ones(1,ns));
%    XZmc=transform_domain(XYpts(:,[1,3]),bl,bu,-s*ones(1,ns),s*ones(1,ns));
XYpts=YZmc;
bl=byz(1,:);
bu=byz(2,:);
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Xb=XY;

 % for xy
%  bl=[min(Xb(:,1))-2,min(Xb(:,2))-2];
%  bu=[max(Xb(:,1))+2,max(Xb(:,2))+2];

%  bl=[min(Xmc(:,1))-2,min(Xmc(:,2))-2,min(Xmc(:,3))-2];
%  bu=[max(Xmc(:,1))+2,max(Xmc(:,2))+2,max(Xmc(:,3))+2];


 s=1;
 XY=transform_domain(Xb,bl,bu,-s*ones(1,ns),s*ones(1,ns));

 
%  
 X=zeros(1,ns,N);
 X(1,:,:)=XY(1:N,:)';
 [y1,M1]=Evol_moments_samples(X,w,1,'raw');
 [y2,M2]=Evol_moments_samples(X,w,2,'raw');
 [y3,M3]=Evol_moments_samples(X,w,3,'raw');
 [y4,M4]=Evol_moments_samples(X,w,4,'raw');
%  [y5,M5]=Evol_moments_samples(X,w,5,'raw');
%  [y6,M6]=Evol_moments_samples(X,w,6,'raw');
 
%  options = statset('MaxIter',5000); 
%  obj = gmdistribution.fit([XY(1:N,1),XY(1:N,2)],20,'Regularize',0.1,'Options',options);
%  
%  [xx,zz]=meshgrid(M1(1)-10*sqrt(M2(1)):0.3:M1(1)+10*sqrt(M2(1)),M1(2)-10*sqrt(M2(3)):0.3:M1(2)+10*sqrt(M2(3)));
%  p=zeros(size(xx));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%  p(i,j)=pdf(obj,[xx(i,j),zz(i,j)]);
%     end
%  end
%  figure
%  plot(XY(1:N,1),XY(1:N,2),'ro')
%  xlabel('x')
%  ylabel('z')
%  hold on
%   contour(xx,zz,p,linspace(0.18,0.001,20))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
y=[zeros(1,ns);y1;y2];
  M=[1;M1';M2'];
   [y,lam,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns))
   save('MaxEnt2_TBP_UTyz','y','lam','bl','bu')
   

%     y=[zeros(1,ns);y1;y2;y3];
%   M=[1;M1';M2';M3'];
%     [y,lam,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns))
%     save('MaxEnt3_TBP_CUT8xz','y','lam','bl','bu')
%     
%     
%     y=[zeros(1,ns);y1;y2;y3;y4];
%   M=[1;M1';M2';M3';M4'];
%   [y,lam,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns))
%   save('MaxEnt4_TBP_CUT8xz','y','lam','bl','bu')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [y,lam,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns))
  
  [Xq,W] = GLeg_pts(29, -s, s);
  W=2^ns*W;

  %% x-y plot
  c=12;
  
   [xx,zz]=meshgrid(linspace(-1*s,1*s,200),linspace(-1*s,1*s,200));
   pent=zeros(size(xx));
 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent(i,j)= pent(i,j)+pdf_MaxEnt([xx(i,j),zz(i,j),Xq(k)],lam,y);
        end
    end
 end
%  c=pdf_MaxEnt(M1,lam,y);
 figure
 plot(XYpts(:,1),XYpts(:,2),'ro')
 xlabel('x')
 ylabel('z')
 hold on
  contour(xx,zz,pent,linspace(0.01,12,20))
   plot_prop_paper
 axis([-s,s,-s,s])
  
  
  figure
   mesh(xx,zz,pent)
 hold on
 xlabel('x')
 ylabel('z')
  plot_prop_paper
axis([-s,s,-s,s])
  
    %% y-z plot
   [xx,zz]=meshgrid(linspace(-1*s,1*s,200),linspace(-1*s,1*s,200));
   pent=zeros(size(xx));
 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([Xq(k),xx(i,j),zz(i,j)],lam,y);
        end
    end
 end
%  c=pdf_MaxEnt(M1,lam,y);
 figure
 plot(XYpts(:,2),XYpts(:,3),'ro')
 xlabel('y')
 ylabel('z')
 hold on
  contour(xx,zz,pent,linspace(0.01,c/20,20))
     plot_prop_paper
 axis([-s,s,-s,s])
  
  
  figure
   mesh(xx,zz,pent)
 hold on
 xlabel('y')
 ylabel('z')
   axis([-s,s,-s,s])
   plot_prop_paper   
    %% x-z plot
   [xx,zz]=meshgrid(linspace(-1*s,1*s,200),linspace(-1*s,1*s,200));
   pent=zeros(size(xx));
 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([xx(i,j),Xq(k),zz(i,j)],lam,y);
        end
    end
 end
%  c=pdf_MaxEnt(M1,lam,y);
 figure
 plot(XYpts(:,1),XYpts(:,3),'ro')
 xlabel('x')
 ylabel('z')
 hold on
  contour(xx,zz,pent,linspace(0.00001,c/5,20))
     plot_prop_paper
 axis([-s,s,-s,s])
  
  
  figure
  mesh(xx,zz,pent)
 hold on
 xlabel('x')
 ylabel('z')
   axis([-s,s,-s,s,0,2*c])
      plot_prop_paper