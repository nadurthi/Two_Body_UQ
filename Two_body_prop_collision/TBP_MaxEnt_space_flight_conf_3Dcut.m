%2 Body problem entropy reconstruction for space flight mnechanics
%conference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   3D reconst CUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
% function MaxEnt2Bodypdfs_time699
% for the 2BP at time step 699 only the pdfs are reconstructed for all 
clear
clc
close all
load('TBPonly699')

%% MC points pdf
NS=[1,2,3];

w=wmc;
X=Data_TBP_at699mc(:,NS);
ns=size(X,2);
Nsamp=size(X,1);


[y,mu1]=Cal_moments_samples(X,w,1,'central');
[y,m2]=Cal_moments_samples(X,w,2,'central');

sqP=sqrtm(inv(moms2cov(y,m2)));
Xt=zeros(size(X));
Xtmc=zeros(size(Data_TBP_at699mc(:,NS)));
for i=1:1:Nsamp
Xt(i,:)=sqP*(X(i,:)'-mu1);
end
for i=1:1:size(Data_TBP_at699mc(:,NS))
Xtmc(i,:)=sqP*(Data_TBP_at699mc(i,NS)'-mu1);
end
 s=8;
 figure
 plot(Xt(:,1),Xt(:,2),'bo')
 
% Xt=Xt(:,[1,2]);
% Xtmc=Xtmc(:,[1,2]);
% NS=[1,2];
% ns=size(Xt,2);


% Xt=X;
% Xtmc=Data_TBP_at699mc(:,NS);

s=1;
[y,mu]=Cal_moments_samples(Xt,w,1,'central');
bl=min(Xt,[],1)-(0/2)*abs(min(Xt,[],1)-mu');
bu=max(Xt,[],1)+(0/2)*abs(max(Xt,[],1)-mu');
[Xt,dt]=transform_domain(Xt,bl,bu,-s*ones(1,ns),s*ones(1,ns));
% transform points used for plotting only
[Xtmc,dtmc]=transform_domain(Xtmc,bl,bu,-s*ones(1,ns),s*ones(1,ns));
figure
plot(Xt(:,1),Xt(:,2),'bo')
figure
plot(Xtmc(:,1),Xtmc(:,2),'bo')


[y1,M1]=Cal_moments_samples(Xt,w,1,'raw');
[y2,M2]=Cal_moments_samples(Xt,w,2,'raw');
[y3,M3]=Cal_moments_samples(Xt,w,3,'raw');
[y4,M4]=Cal_moments_samples(Xt,w,4,'raw');
[y5,M5]=Cal_moments_samples(Xt,w,5,'raw');
[y6,M6]=Cal_moments_samples(Xt,w,6,'raw');
%%%%%%%%%%%%%%% 2nd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=[zeros(1,ns);y1;y2];
/1h ''/1n.0
M=[1;M1;M2];
[Y2,lam2,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),zeros(size(y),1));

%%%%%%%%%%%%%%% 3rd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3];
lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3];
[Y3,lam3,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 4th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4];
lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3;M4];
[Y4,lam4,xl,xu]=MaxEntPdf(y,M,-1.2*s*ones(1,ns),1.2*s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 5th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4;y5];
lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3;M4;M5];
[Y5,lam5,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 6th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4;y5;y6];
lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3;M4;M5;M6];
[Y6,lam6,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

%% computing pdf values for pltoiing
%  [xx,zz]=meshgrid(linspace(-1.5*s,1.5*s,200),linspace(-1.5*s,1.5*s,200));
%    pent2=zeros(size(xx));
%    pent3=zeros(size(xx));
%    pent4=zeros(size(xx));
%    pent5=zeros(size(xx));
%    pent6=zeros(size(xx));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%  pent2(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam2,Y2);
%  pent3(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam3,Y3);
%  pent4(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam4,Y4);
%  pent5(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam5,Y5);
%   pent6(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam6,Y6);
%     end
%  end

 %% computing the marginalised pdf values for plotting xy
 [xx,zz]=meshgrid(linspace(-1.5*s,1.5*s,100),linspace(-1.5*s,1.5*s,100));
 Xrt=zeros(size(xx));
 Yrt=zeros(size(xx));
 
   pent2=zeros(size(xx));
   pent3=zeros(size(xx));
   pent4=zeros(size(xx));
   pent5=zeros(size(xx));
   pent6=zeros(size(xx));

%    [XX,W] = GH_points(0,0.02,35);
% W=W./mvnpdf(XX, 0,0.02);
[XX,W] = GLeg_pts(35*ones(size(xl)), -1.5*s,1.5*s);

 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent2(i,j)=pent2(i,j)+W(k)*pdf_MaxEnt([xx(i,j),zz(i,j),XX(k,:)],lam2,Y2);
 pent3(i,j)=pent3(i,j)+W(k)*pdf_MaxEnt([xx(i,j),zz(i,j),XX(k,:)],lam3,Y3);
 pent4(i,j)=pent4(i,j)+W(k)*pdf_MaxEnt([xx(i,j),zz(i,j),XX(k,:)],lam4,Y4);
 pent5(i,j)=pent5(i,j)+W(k)*pdf_MaxEnt([xx(i,j),zz(i,j),XX(k,:)],lam5,Y5);
 pent6(i,j)=pent6(i,j)+W(k)*pdf_MaxEnt([xx(i,j),zz(i,j),XX(k,:)],lam6,Y6);
         end
    end
    i
 end

 figure
 plot(Xtmc(:,1),Xtmc(:,2),'ro')
  hold on
 contour(xx,zz,pent2,linspace(0.001,max(max(pent2)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
  saveas(gcf,'CUT8699_3D_2mom_xy_1.fig')
 saveas(gcf,'CUT8699_3D_2mom_xy_1.eps') 
  
 figure
 surfl(xx,zz,pent2)
   shading interp
colormap(gray);
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'CUT8699_3D_2mom_xy_2.fig')
 saveas(gcf,'CUT8699_3D_2mom_xy_2.eps') 
 
 
  figure
 plot(Xtmc(:,1),Xtmc(:,2),'ro')
 hold on
 contour(xx,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
    saveas(gcf,'CUT8699_3D_3mom_xy_1.fig')
 saveas(gcf,'CUT8699_3D_3mom_xy_1.eps') 
  
 figure
  surfl(xx,zz,pent3)
    shading interp
colormap(gray);
   plot_prop_paper
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
 saveas(gcf,'CUT8699_3D_3mom_xy_2.fig')
 saveas(gcf,'CUT8699_3D_3mom_xy_2.eps') 
 
  figure
  plot(Xtmc(:,1),Xtmc(:,2),'ro')
   hold on
 contour(xx,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'CUT8699_3D_4mom_xy_1.fig')
 saveas(gcf,'CUT8699_3D_4mom_xy_1.eps') 
 
 figure
  surfl(xx,zz,pent4)
    shading interp
colormap(gray);
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
 saveas(gcf,'CUT8699_3D_4mom_xy_2.fig')
 saveas(gcf,'CUT8699_3D_4mom_xy_2.eps') 
 
 
  figure
  plot(Xtmc(:,1),Xtmc(:,2),'ro')
   hold on
 contour(xx,zz,pent5,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
     xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
  saveas(gcf,'CUT8699_3D_5mom_xy_1.fig')
 saveas(gcf,'CUT8699_3D_5mom_xy_1.eps') 
 
 figure
  surfl(xx,zz,pent5)
    shading interp
colormap(gray);
      xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
   saveas(gcf,'CUT8699_3D_5mom_xy_2.fig')
 saveas(gcf,'CUT8699_3D_5mom_xy_2.eps')
 
 
  figure
  plot(Xtmc(:,1),Xtmc(:,2),'ro')
   hold on
 contour(xx,zz,pent6,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
    xlabel('x')
   ylabel('y')
   legend('MC samples','PME-pdf contours')
  plot_prop_paper
   axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
      saveas(gcf,'CUT8699_3D_6mom_xy_1.fig')
 saveas(gcf,'CUT8699_3D_6mom_xy_1.eps')
 
 figure
  surfl(xx,zz,pent6)
  shading interp
colormap(gray);
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
       saveas(gcf,'CUT8699_3D_6mom_xy_2.fig')
 saveas(gcf,'CUT8699_3D_6mom_xy_2.eps')
 %%%%%%%%%%%%%%%%%%%%%%%%%
%  yz - marginalization
 %% computing the marginalised pdf values for plotting xy
 [yy,zz]=meshgrid(linspace(-1.5*s,1.5*s,100),linspace(-1.5*s,1.5*s,100));
 
   pent2=zeros(size(xx));
   pent3=zeros(size(xx));
   pent4=zeros(size(xx));
   pent5=zeros(size(xx));
   pent6=zeros(size(xx));

%    [XX,W] = GH_points(0,0.02,45);
% W=W./mvnpdf(XX, 0,0.02);
% [XX,W] = GLeg_pts(29*ones(size(xl)), -2*s,2*s);

 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent2(i,j)=pent2(i,j)+W(k)*pdf_MaxEnt([XX(k,:),yy(i,j),zz(i,j)],lam2,Y2);
 pent3(i,j)=pent3(i,j)+W(k)*pdf_MaxEnt([XX(k,:),yy(i,j),zz(i,j)],lam3,Y3);
 pent4(i,j)=pent4(i,j)+W(k)*pdf_MaxEnt([XX(k,:),yy(i,j),zz(i,j)],lam4,Y4);
 pent5(i,j)=pent5(i,j)+W(k)*pdf_MaxEnt([XX(k,:),yy(i,j),zz(i,j)],lam5,Y5);
 pent6(i,j)=pent6(i,j)+W(k)*pdf_MaxEnt([XX(k,:),yy(i,j),zz(i,j)],lam6,Y6);
         end
    end
    i
 end
 %% plotting

 figure
 plot(Xtmc(:,2),Xtmc(:,3),'ro')
  hold on
 contour(yy,zz,pent2,linspace(0.001,max(max(pent2)/3),15),'linewidth',2)
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'CUT8699_3D_2mom_yz_1.fig')
 saveas(gcf,'CUT8699_3D_2mom_yz_1.eps')
  
 figure
 surfl(yy,zz,pent2)
   shading interp
colormap(gray);
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
    saveas(gcf,'CUT8699_3D_2mom_yz_2.fig')
 saveas(gcf,'CUT8699_3D_2mom_yz_2.eps')
 %%
  figure
 plot(Xtmc(:,2),Xtmc(:,3),'ro')
 hold on
 contour(yy,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'CUT8699_3D_3mom_yz_1.fig')
 saveas(gcf,'CUT8699_3D_3mom_yz_1.eps') 
  
 figure
  surfl(yy,zz,pent3)
    shading interp
colormap(gray);
   plot_prop_paper
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'CUT8699_3D_3mom_yz_2.fig')
 saveas(gcf,'CUT8699_3D_3mom_yz_2.eps') 
 %%
 
  figure
  plot(Xtmc(:,2),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
   saveas(gcf,'CUT8699_3D_4mom_yz_1.fig')
 saveas(gcf,'CUT8699_3D_4mom_yz_1.eps') 
  
 figure
  surfl(yy,zz,pent4)
    shading interp
colormap(gray);
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
   saveas(gcf,'CUT8699_3D_4mom_yz_2.fig')
 saveas(gcf,'CUT8699_3D_4mom_yz_2.eps') 

%%
  figure
  plot(Xtmc(:,2),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent5,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
     xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'CUT8699_3D_5mom_yz_1.fig')
 saveas(gcf,'CUT8699_3D_5mom_yz_1.eps') 

 
 figure
  surfl(yy,zz,pent5)
    shading interp
colormap(gray);
      xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'CUT8699_3D_5mom_yz_2.fig')
 saveas(gcf,'CUT8699_3D_5mom_yz_2.eps') 
 
 %% 
  figure
  plot(Xtmc(:,2),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent6,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
    xlabel('y')
   ylabel('z')
   legend('MC samples','PME-pdf contours')
  plot_prop_paper
   axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
    saveas(gcf,'CUT8699_3D_6mom_yz_1.fig')
 saveas(gcf,'CUT8699_3D_6mom_yz_1.eps') 
 
 figure
  surfl(yy,zz,pent6)
  shading interp
colormap(gray);
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
     saveas(gcf,'CUT8699_3D_6mom_yz_2.fig')
 saveas(gcf,'CUT8699_3D_6mom_yz_2.eps') 
 
  %%%%%%%%%%%%%%%%%%%%%%%%%
%  xz - marginalization
 %% computing the marginalised pdf values for plotting 
 [yy,zz]=meshgrid(linspace(-1.5*s,1.5*s,100),linspace(-1.5*s,1.5*s,100));
 
   pent2=zeros(size(xx));
   pent3=zeros(size(xx));
   pent4=zeros(size(xx));
   pent5=zeros(size(xx));
   pent6=zeros(size(xx));

%    [XX,W] = GH_points(0,0.02,45);
% W=W./mvnpdf(XX, 0,0.02);
% [XX,W] = GLeg_pts(29*ones(size(xl)), -2*s,2*s);

 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent2(i,j)=pent2(i,j)+W(k)*pdf_MaxEnt([yy(i,j),XX(k,:),zz(i,j)],lam2,Y2);
 pent3(i,j)=pent3(i,j)+W(k)*pdf_MaxEnt([yy(i,j),XX(k,:),zz(i,j)],lam3,Y3);
 pent4(i,j)=pent4(i,j)+W(k)*pdf_MaxEnt([yy(i,j),XX(k,:),zz(i,j)],lam4,Y4);
 pent5(i,j)=pent5(i,j)+W(k)*pdf_MaxEnt([yy(i,j),XX(k,:),zz(i,j)],lam5,Y5);
 pent6(i,j)=pent6(i,j)+W(k)*pdf_MaxEnt([yy(i,j),XX(k,:),zz(i,j)],lam6,Y6);
         end
    end
    i
 end
 %% plotting

 figure
 plot(Xtmc(:,1),Xtmc(:,3),'ro')
  hold on
 contour(yy,zz,pent2,linspace(0.001,max(max(pent2)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'CUT8699_3D_2mom_xz_1.fig')
 saveas(gcf,'CUT8699_3D_2mom_xz_1.eps')
  
 figure
 surfl(yy,zz,pent2)
   shading interp
colormap(gray);
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'CUT8699_3D_2mom_xz_2.fig')
 saveas(gcf,'CUT8699_3D_2mom_xz_2.eps')
 %%
  figure
 plot(Xtmc(:,1),Xtmc(:,3),'ro')
 hold on
 contour(yy,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
   saveas(gcf,'CUT8699_3D_3mom_xz_1.fig')
 saveas(gcf,'CUT8699_3D_3mom_xz_1.eps')
  
 figure
  surfl(yy,zz,pent3)
    shading interp
colormap(gray);
   plot_prop_paper
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
    saveas(gcf,'CUT8699_3D_3mom_xz_2.fig')
 saveas(gcf,'CUT8699_3D_3mom_xz_2.eps')
 %%
 
  figure
  plot(Xtmc(:,1),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
     saveas(gcf,'CUT8699_3D_4mom_xz_1.fig')
 saveas(gcf,'CUT8699_3D_4mom_xz_1.eps')
 
 
 figure
  surfl(yy,zz,pent4)
    shading interp
colormap(gray);
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
    saveas(gcf,'CUT8699_3D_4mom_xz_2.fig')
 saveas(gcf,'CUT8699_3D_4mom_xz_2.eps')
%%
  figure
  plot(Xtmc(:,1),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent5,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
     xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
     saveas(gcf,'CUT8699_3D_5mom_xz_1.fig')
 saveas(gcf,'CUT8699_3D_5mom_xz_1.eps')
 
 
 figure
 mesh(yy,zz,pent5)
  surfl(yy,zz,pent5)
    shading interp
colormap(gray);
      xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
      saveas(gcf,'CUT8699_3D_5mom_xz_2.fig')
 saveas(gcf,'CUT8699_3D_5mom_xz_2.eps')
 
 %% 
  figure
  plot(Xtmc(:,1),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent6,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
    xlabel('x')
   ylabel('z')
   legend('MC samples','PME-pdf contours')
  plot_prop_paper
   axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
        saveas(gcf,'CUT8699_3D_6mom_xz_1.fig')
 saveas(gcf,'CUT8699_3D_6mom_xz_1.eps')
 
 
 figure
  surfl(yy,zz,pent6)
  shading interp
colormap(gray);
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
 saveas(gcf,'CUT8699_3D_6mom_xz_2.fig')
 saveas(gcf,'CUT8699_3D_6mom_xz_2.eps')
 close all
 %%%%%%%%%%%%%%%%%%%%%%%%%
%  %% now plotting the transformed pdf to the original config
% [y,m2]=Cal_moments_samples(X,w,2,'central');
% sqP=sqrtm(inv(moms2cov(y,m2)));
% P=sqrtm(moms2cov(y,m2));
% b=max(diag(P));
% [xx,zz]=meshgrid(linspace(mu1(1)-b,mu1(1)+b,100),linspace(mu1(2)-b,mu1(2)+b,100));
%   
%    pent2t=zeros(size(xx));
%    pent3t=zeros(size(xx));
%    pent4t=zeros(size(xx));
%    pent5t=zeros(size(xx));
%    pent6t=zeros(size(xx));
% 
% %    [XX,W] = GH_points(0,1,45);
% % W=W./mvnpdf(XX, 0, 1);
% [XX,W] = GLeg_pts(29*ones(size(xl)),mu1(3)-b,mu1(3)+b);
% C=det(inv(sqP)*diag(1./dt));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         for k=1:1:length(W)
%             XY=sqP*([xx(i,j),zz(i,j),XX(k,:)]'-mu1);
%             [XY,dt]=transform_domain(XY',bl,bu,-s*ones(1,ns),s*ones(1,ns));
%  pent2t(i,j)=pent2t(i,j)+W(k)*pdf_MaxEnt_extd(XY,xl,xu,lam2,Y2)/C;
%  pent3t(i,j)=pent3t(i,j)+W(k)*pdf_MaxEnt_extd(XY,xl,xu,lam3,Y3)/C;
%  pent4t(i,j)=pent4t(i,j)+W(k)*pdf_MaxEnt_extd(XY,xl,xu,lam4,Y4)/C;
%  pent5t(i,j)=pent5t(i,j)+W(k)*pdf_MaxEnt_extd(XY,xl,xu,lam5,Y5)/C;
%  pent6t(i,j)=pent6t(i,j)+W(k)*pdf_MaxEnt_extd(XY,xl,xu,lam6,Y6)/C;
%          end
%     end
%     i
%  end
%  
%  %% 
%  figure
%  plot(X(:,1),X(:,2),'ro')
%   hold on
% %   plot(Xrt(:,1),Yrt(:,1),'b+')
%  contour(xx,zz,pent2t,15,'linewidth',2)
%    xlabel('x')
%    ylabel('y')
%    plot_prop_paper
% %  axis([-2*s,2*s,-1.5*s,1.5*s])
%   for i=1:1:size(xx,1)
%       for j=1:1:size(xx,2)
%           plot(xx(i,j),zz(i,j),'b+')
%       end
%   end
%   
%  figure
%   plot(X(:,1),X(:,2),'ro')
%   hold on
%  surfl(xx,zz,pent4t)
%    shading interp
% colormap(gray);
%    xlabel('x')
%    ylabel('y')
%    plot_prop_paper
%  axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
%  view([44,44])