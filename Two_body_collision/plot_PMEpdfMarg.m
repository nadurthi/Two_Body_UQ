function plot_PMEpdfMarg(YY,R)
% S=[MD.M2(1),MD.M2(3),MD.M2(2);MD.M2(3),MD.M2(6),MD.M2(4);MD.M2(2),MD.M2(4),MD.M2(5)];
% S=S-muu'*muu;
S=diag([ 1.0791e-06, 0.0001150, 0.004]);
muu=YY.muu';
%% x-y
[xx,yy]=meshgrid(linspace(muu(1)-3*sqrt(S(1,1)),muu(1)+3*sqrt(S(1,1)),100),linspace(muu(2)-3*sqrt(S(2,2)),muu(2)+3*sqrt(S(2,2)),100));
pent=zeros(size(xx));

[XX,W] = GH_points(muu(3),S(3,3),35);
W=W./mvnpdf(XX, muu(3),S(3,3));

iP=YY.iP;
 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
        sp=(iP*([xx(i,j),yy(i,j),XX(k,:)]-muu)')';    
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],YY.lam,YY.Y)*det(iP);
         end
    end
 end
figure(1)
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 hold on
 contour(xx,yy,pent,15)
 xlabel('x')
 ylabel('y')
%  plot_prop_paper
%  saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'fig')
  
 figure(2)
 mesh(xx,yy,pent)
 hold on
%  plot([R,R,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k--','linewidth',3)
 xlabel('x')
 ylabel('y')
%   plot_prop_paper
%  saveas(gcf,strcat('MISSxy2_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSxy2_',num2str(kkk)),'fig')
  

  %% y-z 

  keyboard
  S=diag([2.0791e-06, 0.0001150, 0.002]);
    S=diag([1, 4, 4]);
% [yy,zz]=meshgrid(linspace(muu(2)-3*sqrt(S(2,2)),muu(2)+3*sqrt(S(2,2)),100),linspace(muu(3)-3*sqrt(S(3,3)),muu(3)+3*sqrt(S(3,3)),100));
% pent=zeros(size(yy));
% [XX,W] = GH_points(muu(1),S(1,1),35);
% W=W./mvnpdf(XX, muu(1),S(1,1));

[yy,zz]=meshgrid(linspace(0-3*sqrt(S(2,2)),0+3*sqrt(S(2,2)),100),linspace(0-3*sqrt(S(3,3)),0+3*sqrt(S(3,3)),100));
pent=zeros(size(yy));

[XX,W] = GH_points(0,S(1,1),35);
W=W./mvnpdf(XX,0,S(1,1));


 for i=1:1:size(yy,1)
    for j=1:1:size(yy,2)
        for k=1:1:length(W)
%         sp=(iPs*([XX(k,:),yy(i,j),zz(i,j)]-muu)')';    
%         pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],YY.lam,YY.Y)*det(iP);
          sp=[XX(k,:),yy(i,j),zz(i,j)];
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],YY.lam,YY.Y);
         end
    end
 end
figure(3)
pp=sort(min(pent));
 contour(yy,zz,pent,linspace(pp(1)*1,max(max(pent)),25))
 hold on
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
  xlabel('y')
 ylabel('z')
%   plot_prop_paper
%   saveas(gcf,strcat('MISSyz1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSyz1_',num2str(kkk)),'fig')
 
 figure(4)
mesh(yy,zz,pent)
 hold on
%   plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 xlabel('y')
 ylabel('z')
%   plot_prop_paper
%  saveas(gcf,strcat('MISSyz2_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSyz2_',num2str(kkk)),'fig')
 
 %% x-z 

[xx,zz]=meshgrid(linspace(muu(1)-3*sqrt(S(1,1)),muu(1)+3*sqrt(S(1,1)),100),linspace(muu(3)-3*sqrt(S(3,3)),muu(3)+3*sqrt(S(3,3)),100));
pent=zeros(size(xx));


[XX,W] = GH_points(muu(2),S(2,2),35);
W=W./mvnpdf(XX, muu(2),S(2,2));


 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
        sp=(iP*([xx(i,j),XX(k,:),zz(i,j)]-muu)')';    
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],YY.lam,YY.Y)*det(iP);
         end
    end
 end
figure(5)
pp=sort(min(pent));
 contour(xx,zz,pent,linspace(pp(2)*1,max(max(pent)),25))
 hold on
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 xlabel('x')
 ylabel('z')
%   plot_prop_paper
%  saveas(gcf,strcat('MISSxz1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSxz1_',num2str(kkk)),'fig')
 
 figure(6)
 mesh(xx,zz,pent)
  hold on
%   plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 xlabel('x')
 ylabel('z')
%   plot_prop_paper
%    saveas(gcf,strcat('MISSxz2_',dtstr),'pdf')
%  saveas(gcf,strcat('MISSxz2_',dtstr),'fig')
%  
%  pause(1)
%  close all
