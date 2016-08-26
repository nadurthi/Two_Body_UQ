%2 Body problem entropy reconstruction
% function MaxEnt2Bodypdfs_polar
% clear
% close all
% clc
% 
%  distcomp.feature( 'LocalUseMpiexec', false );
% matlabpool open 8
w=wmc;
X=XYmc;
ns=size(X,2);
Nsamp=size(X,1);


% [y,mu1]=Cal_moments_samples(X,w,1,'central');
% [y,m2]=Cal_moments_samples(X,w,2,'central');
% 
% sqP=sqrtm(inv(moms2cov(y,m2)));
% Xt=zeros(size(X));
% Xtmc=zeros(size(Data_TBP_at699mc(:,NS)));
% for i=1:1:Nsamp
% Xt(i,:)=sqP*(X(i,:)'-mu1);
% end
% for i=1:1:size(Data_TBP_at699mc(:,NS))
% Xtmc(i,:)=sqP*(Data_TBP_at699mc(i,NS)'-mu1);
% end
%  s=8;
%  figure
%  plot(Xt(:,1),Xt(:,2),'bo')
 
% Xt=Xt(:,[1,2]);
% Xtmc=Xtmc(:,[1,2]);
% NS=[1,2];
% ns=size(Xt,2);


Xt=X;
Xtmc=XYmc;

s=1;
[y,mu]=Cal_moments_samples(Xt,w,1,'central');
bl=1*min(Xt,[],1)-(1/2)*abs(min(Xt,[],1)-mu');
bu=1*max(Xt,[],1)+(1/2)*abs(max(Xt,[],1)-mu');
[Xt,dt]=transform_domain(Xt,bl,bu,-s*ones(1,ns),s*ones(1,ns));
% transform points used for plotting only
[Xtmc,dtmc]=transform_domain(Xtmc,bl,bu,-s*ones(1,ns),s*ones(1,ns));
figure
% plot(Xt(:,1),Xt(:,2),'bo')
axis([-1,1,-1,1])


[y1,M1]=Cal_moments_samples(Xt,w,1,'raw');
[y2,M2]=Cal_moments_samples(Xt,w,2,'raw');
[y3,M3]=Cal_moments_samples(Xt,w,3,'raw');
[y4,M4]=Cal_moments_samples(Xt,w,4,'raw');
% [y5,M5]=Cal_moments_samples(Xt,w,5,'raw');
% [y6,M6]=Cal_moments_samples(Xt,w,6,'raw');
%%%%%%%%%%%%%%% 2nd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y=[zeros(1,ns);y1;y2];

% M=[1;M1;M2];
% [Y2,lam2,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),zeros(size(y),1));

%%%%%%%%%%%%%%% 3rd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y=[zeros(1,ns);y1;y2;y3];
% lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
% M=[1;M1;M2;M3];
% [Y3,lam3,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 4th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,ns);y1;y2;y3;y4];
lam0=zeros(size(y,1),1);
M=[1;M1;M2;M3;M4];
[Y4,lam4,xl,xu]=MaxEntPdf(y,M,-s*ones(1,ns),s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 5th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y=[zeros(1,ns);y1;y2;y3;y4;y5];
% lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
% M=[1;M1;M2;M3;M4;M5];
% [Y5,lam5,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 6th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y=[zeros(1,ns);y1;y2;y3;y4;y5;y6];
% lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
% M=[1;M1;M2;M3;M4;M5;M6];
% [Y6,lam6,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

% %% computing pdf values for pltoiing
 [xx,zz]=meshgrid(linspace(-s,1*s,200),linspace(-1*s,1*s,200));
   pent2=zeros(size(xx));
   pent3=zeros(size(xx));
   pent4=zeros(size(xx));
%    pent5=zeros(size(xx));
%    pent6=zeros(size(xx));
 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
%  pent2(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam2,Y2);
%  pent3(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam3,Y3);
 pent4(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam4,Y4);
%  pent5(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam5,Y5);
%   pent6(i,j)=pdf_MaxEnt([xx(i,j),zz(i,j)],lam6,Y6);
    end
 end
%  figure
%  plot(Xtmc(:,1),Xtmc(:,2),'ro')
%   hold on
%  contour(xx,zz,pent2,linspace(0.001,max(max(pent2)/3),15),'linewidth',2)
%   if NS(1)==1
%       xlabel('x')
%   elseif NS(1)==2
%       xlabel('y')
%   else
%       xlabel('z')
%   end
%     if NS(2)==1
%       ylabel('x')
%   elseif NS(2)==2
%       ylabel('y')
%   else
%       ylabel('z')
%     end
%   
%    plot_prop_paper
%  axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
%   saveas(gcf,strcat(name,'_2mom_1.fig'))
%   saveas(gcf,strcat(name,'_2mom_1.eps'))
%   
%  figure
%  surfl(xx,zz,pent2)
%    shading interp
% colormap(gray);
%   if NS(1)==1
%       xlabel('x')
%   elseif NS(1)==2
%       xlabel('y')
%   else
%       xlabel('z')
%   end
%     if NS(2)==1
%       ylabel('x')
%   elseif NS(2)==2
%       ylabel('y')
%   else
%       ylabel('z')
%     end
%    plot_prop_paper
%  axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
%  view([75,56])
%   saveas(gcf,strcat(name,'_2mom_2.fig'))
%   saveas(gcf,strcat(name,'_2mom_2.eps'))
%   
%   %
%   figure
%  plot(Xtmc(:,1),Xtmc(:,2),'ro')
%   hold on
%  contour(xx,zz,pent3,linspace(0.001,max(max(pent3)/3),15),'linewidth',2)
%   if NS(1)==1
%       xlabel('x')
%   elseif NS(1)==2
%       xlabel('y')
%   else
%       xlabel('z')
%   end
%     if NS(2)==1
%       ylabel('x')
%   elseif NS(2)==2
%       ylabel('y')
%   else
%       ylabel('z')
%     end
%    plot_prop_paper
%  axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
%   saveas(gcf,strcat(name,'_3mom_1.fig'))
%   saveas(gcf,strcat(name,'_3mom_1.eps'))
%   
%  figure
%  surfl(xx,zz,pent3)
%    shading interp
% colormap(gray);
%   if NS(1)==1
%       xlabel('x')
%   elseif NS(1)==2
%       xlabel('y')
%   else
%       xlabel('z')
%   end
%     if NS(2)==1
%       ylabel('x')
%   elseif NS(2)==2
%       ylabel('y')
%   else
%       ylabel('z')
%     end
%    plot_prop_paper
%  axis([-1.5*s,1.5*s,-1.5*s,1.5*s])
% view([75,56])
%   saveas(gcf,strcat(name,'_3mom_2.fig'))
%   saveas(gcf,strcat(name,'_3mom_2.eps'))
%   
%   
  %
  figure
 plot(Xtmc(:,1),Xtmc(:,2),'ro')
  hold on
 contour(xx,zz,pent4,linspace(0.001,max(max(pent4)/3),15),'linewidth',2)
%   if NS(1)==1
%       xlabel('x')
%   elseif NS(1)==2
%       xlabel('y')
%   else
%       xlabel('z')
%   end
%     if NS(2)==1
%       ylabel('x')
%   elseif NS(2)==2
%       ylabel('y')
%   else
%       ylabel('z')
%     end
   plot_prop_paper
 axis([-1*s,1*s,-1*s,1*s])
%   saveas(gcf,strcat(name,'_4mom_1.fig'))
%   saveas(gcf,strcat(name,'_4mom_1.eps'))
  