%2 Body problem entropy reconstruction for space flight mnechanics
%conference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   3D reconst %%%%%%%%%%%%%%%%%%%%%%%%%%%
% function MaxEnt2Bodypdfs_time699
% for the 2BP at time step 699 only the pdfs are reconstructed for all 
function TBP_maxent_abramov
clear
clc
close all
load('TBPonly699')
global pic_cnt
pic_cnt=1;
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



[y1,M1]=Cal_moments_samples(Xt,w,1,'raw');
[y2,M2]=Cal_moments_samples(Xt,w,2,'raw');
[y3,M3]=Cal_moments_samples(Xt,w,3,'raw');
[y4,M4]=Cal_moments_samples(Xt,w,4,'raw');
[y5,M5]=Cal_moments_samples(Xt,w,5,'raw');
[y6,M6]=Cal_moments_samples(Xt,w,6,'raw');

xl=-1.5*s*ones(1,ns);
xu=1.5*s*ones(1,ns);
method.xl=xl;
method.xu=xu;
method.name='gl';
method.glN=25;
% DD=@(y)magic(size(y,1))'*magic(size(y,1));
DD=@(y)eye(size(y,1));
%%%%%%%%%%%%%%% 2nd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=[zeros(1,ns);y1;y2];
[V,D]=eig(DD(y));
method.A=V;

M=[1;M1;M2];
[Y2,lam2]=MaxEntPdf_abramov(y,M,zeros(size(y),1),method);

%%%%%%%%%%%%%%% 3rd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3];
[V,D]=eig(DD(y));
method.A=V;

lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3];
[Y3,lam3]=MaxEntPdf_abramov(y,M,lam0,method);

%%%%%%%%%%%%%%% 4th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4];
[V,D]=eig(DD(y));
method.A=V;

lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3;M4];
[Y4,lam4]=MaxEntPdf_abramov(y,M,lam0,method);

%%%%%%%%%%%%%%% 5th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4;y5];
[V,D]=eig(DD(y));
method.A=V;

lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3;M4;M5];
[Y5,lam5]=MaxEntPdf_abramov(y,M,lam0,method);

%%%%%%%%%%%%%%% 6th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4;y5;y6];
[V,D]=eig(DD(y));
method.A=V;

lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
M=[1;M1;M2;M3;M4;M5;M6];
[Y6,lam6]=MaxEntPdf_abramov(y,M,lam0,method);

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
pause(0.5)
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
[XX,W] = GLeg_pts(29*ones(size(xl)), -1.5*s,1.5*s);

 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
 pent2(i,j)=pent2(i,j)+W(k)*pdf_MaxEnt_extd([xx(i,j),zz(i,j),XX(k,:)],lam2,Y2,xl,xu,A);
 pent3(i,j)=pent3(i,j)+W(k)*pdf_MaxEnt_extd([xx(i,j),zz(i,j),XX(k,:)],lam3,Y3,xl,xu,A);
 pent4(i,j)=pent4(i,j)+W(k)*pdf_MaxEnt_extd([xx(i,j),zz(i,j),XX(k,:)],lam4,Y4,xl,xu,A);
 pent5(i,j)=pent5(i,j)+W(k)*pdf_MaxEnt_extd([xx(i,j),zz(i,j),XX(k,:)],lam5,Y5,xl,xu,A);
 pent6(i,j)=pent6(i,j)+W(k)*pdf_MaxEnt_extd([xx(i,j),zz(i,j),XX(k,:)],lam6,Y6,xl,xu,A);
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
  saveas(gcf,'MC699_3D_2mom_xy_1.fig')
 saveas(gcf,'MC699_3D_2mom_xy_1.eps') 
  
 figure
 surfl(xx,zz,pent2)
   shading interp
colormap(gray);
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'MC699_3D_2mom_xy_2.fig')
 saveas(gcf,'MC699_3D_2mom_xy_2.eps') 
 
 
  figure
 plot(Xtmc(:,1),Xtmc(:,2),'ro')
 hold on
 contour(xx,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
    saveas(gcf,'MC699_3D_3mom_xy_1.fig')
 saveas(gcf,'MC699_3D_3mom_xy_1.eps') 
  
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
 saveas(gcf,'MC699_3D_3mom_xy_2.fig')
 saveas(gcf,'MC699_3D_3mom_xy_2.eps') 
 
  figure
  plot(Xtmc(:,1),Xtmc(:,2),'ro')
   hold on
 contour(xx,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'MC699_3D_4mom_xy_1.fig')
 saveas(gcf,'MC699_3D_4mom_xy_1.eps') 
 
 figure
  surfl(xx,zz,pent4)
    shading interp
colormap(gray);
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
 saveas(gcf,'MC699_3D_4mom_xy_2.fig')
 saveas(gcf,'MC699_3D_4mom_xy_2.eps') 
 
 
  figure
  plot(Xtmc(:,1),Xtmc(:,2),'ro')
   hold on
 contour(xx,zz,pent5,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
     xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
  saveas(gcf,'MC699_3D_5mom_xy_1.fig')
 saveas(gcf,'MC699_3D_5mom_xy_1.eps') 
 
 figure
  surfl(xx,zz,pent5)
    shading interp
colormap(gray);
      xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
   saveas(gcf,'MC699_3D_5mom_xy_2.fig')
 saveas(gcf,'MC699_3D_5mom_xy_2.eps')
 
 
  figure
  plot(Xtmc(:,1),Xtmc(:,2),'ro')
   hold on
 contour(xx,zz,pent6,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
    xlabel('x')
   ylabel('y')
   legend('MC samples','PME-pdf contours')
  plot_prop_paper
   axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
      saveas(gcf,'MC699_3D_6mom_xy_1.fig')
 saveas(gcf,'MC699_3D_6mom_xy_1.eps')
 
 figure
  surfl(xx,zz,pent6)
  shading interp
colormap(gray);
   xlabel('x')
   ylabel('y')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
       saveas(gcf,'MC699_3D_6mom_xy_2.fig')
 saveas(gcf,'MC699_3D_6mom_xy_2.eps')
 pause(0.5)
 %%%%%%%%%%%%%%%%%%%%%%%%%
%  yz - marginalization
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
 saveas(gcf,'MC699_3D_2mom_yz_1.fig')
 saveas(gcf,'MC699_3D_2mom_yz_1.eps')
  
 figure
 surfl(yy,zz,pent2)
   shading interp
colormap(gray);
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
    saveas(gcf,'MC699_3D_2mom_yz_2.fig')
 saveas(gcf,'MC699_3D_2mom_yz_2.eps')
 %%
  figure
 plot(Xtmc(:,2),Xtmc(:,3),'ro')
 hold on
 contour(yy,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'MC699_3D_3mom_yz_1.fig')
 saveas(gcf,'MC699_3D_3mom_yz_1.eps') 
  
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
  saveas(gcf,'MC699_3D_3mom_yz_2.fig')
 saveas(gcf,'MC699_3D_3mom_yz_2.eps') 
 %%
 
  figure
  plot(Xtmc(:,2),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
   saveas(gcf,'MC699_3D_4mom_yz_1.fig')
 saveas(gcf,'MC699_3D_4mom_yz_1.eps') 
  
 figure
  surfl(yy,zz,pent4)
    shading interp
colormap(gray);
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
   saveas(gcf,'MC699_3D_4mom_yz_2.fig')
 saveas(gcf,'MC699_3D_4mom_yz_2.eps') 

%%
  figure
  plot(Xtmc(:,2),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent5,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
     xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 saveas(gcf,'MC699_3D_5mom_yz_1.fig')
 saveas(gcf,'MC699_3D_5mom_yz_1.eps') 

 
 figure
  surfl(yy,zz,pent5)
    shading interp
colormap(gray);
      xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'MC699_3D_5mom_yz_2.fig')
 saveas(gcf,'MC699_3D_5mom_yz_2.eps') 
 
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
    saveas(gcf,'MC699_3D_6mom_yz_1.fig')
 saveas(gcf,'MC699_3D_6mom_yz_1.eps') 
 
 figure
  surfl(yy,zz,pent6)
  shading interp
colormap(gray);
   xlabel('y')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
     saveas(gcf,'MC699_3D_6mom_yz_2.fig')
 saveas(gcf,'MC699_3D_6mom_yz_2.eps') 
 pause(0.5)
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
 saveas(gcf,'MC699_3D_2mom_xz_1.fig')
 saveas(gcf,'MC699_3D_2mom_xz_1.eps')
  
 figure
 surfl(yy,zz,pent2)
   shading interp
colormap(gray);
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
  saveas(gcf,'MC699_3D_2mom_xz_2.fig')
 saveas(gcf,'MC699_3D_2mom_xz_2.eps')
 %%
  figure
 plot(Xtmc(:,1),Xtmc(:,3),'ro')
 hold on
 contour(yy,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
   saveas(gcf,'MC699_3D_3mom_xz_1.fig')
 saveas(gcf,'MC699_3D_3mom_xz_1.eps')
  
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
    saveas(gcf,'MC699_3D_3mom_xz_2.fig')
 saveas(gcf,'MC699_3D_3mom_xz_2.eps')
 %%
 
  figure
  plot(Xtmc(:,1),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
     saveas(gcf,'MC699_3D_4mom_xz_1.fig')
 saveas(gcf,'MC699_3D_4mom_xz_1.eps')
 
 
 figure
  surfl(yy,zz,pent4)
    shading interp
colormap(gray);
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
    saveas(gcf,'MC699_3D_4mom_xz_2.fig')
 saveas(gcf,'MC699_3D_4mom_xz_2.eps')
%%
  figure
  plot(Xtmc(:,1),Xtmc(:,3),'ro')
   hold on
 contour(yy,zz,pent5,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
     xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
     saveas(gcf,'MC699_3D_5mom_xz_1.fig')
 saveas(gcf,'MC699_3D_5mom_xz_1.eps')
 
 
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
      saveas(gcf,'MC699_3D_5mom_xz_2.fig')
 saveas(gcf,'MC699_3D_5mom_xz_2.eps')
 
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
        saveas(gcf,'MC699_3D_6mom_xz_1.fig')
 saveas(gcf,'MC699_3D_6mom_xz_1.eps')
 
 
 figure
  surfl(yy,zz,pent6)
  shading interp
colormap(gray);
   xlabel('x')
   ylabel('z')
   plot_prop_paper
 axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
 view([44,44])
 saveas(gcf,'MC699_3D_6mom_xz_2.fig')
 saveas(gcf,'MC699_3D_6mom_xz_2.eps')
 pause(0.5)
end
%%
function [y,lam]=MaxEntPdf_abramov(y,M,lam0,method)
ns=size(y,2);
nm=size(y,1);
options=optimset('disp','iter','TolCon',1e-7,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-5,'TolX',1e-5,'Jacobian','on');
lam=fsolve(@(lam)maxentFSOLVE_abramov(lam,y,M,method),lam0,options);
end
%%
function [ceq,jac]=maxentFSOLVE_abramov(lam,y,M,method)
ns=size(y,2);
A=method.A;
xl=method.xl;
xu=method.xu;


if strcmp(method.name,'gh')==1
[X,W] = GH_points(method.ghmu,method.ghP,method.ghN);
W=W./mvnpdf(X, zeros(size(y,2),1)', 0.02*eye(size(y,2)));
end
if strcmp(method.name,'gl')==1
[X,W] = GLeg_pts(method.glN*ones(size(method.xl)), xl,xu);
W=2^ns*W;
end

nm=size(y,1);
nq=length(W);
ceq=zeros(nm,nq);
jac=zeros(nm^2,nq);
for i=1:length(W)
    ceq(:,i)=W(i)*A*prod(repmat(X(i,:),nm,1).^y,2)*pdf_MaxEnt_extd(X(i,:),lam,y,xl,xu,A);
    jac(:,i)=reshape(W(i)*A*(prod(repmat(X(i,:),nm,1).^y,2)*prod(repmat(X(i,:),nm,1).^y,2)')*A'*pdf_MaxEnt_extd(X(i,:),lam,y,xl,xu,A),nm^2,1);
end
ceq=sum(ceq,2)-A*M;
jac=reshape(sum(jac,2),nm,nm);
end