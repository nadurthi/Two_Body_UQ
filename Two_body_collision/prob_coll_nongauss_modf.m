function prob_coll=prob_coll_nongauss_modf(MD,Nmoms,R,kkk)
dim=length(MD.M1);
tic
y1=[     1     0     0
        0     1     0
        0     0     1];
% 
y2=  [2     0     0
     1     0     1
     1     1     0
     0     1     1
     0     0     2
     0     2     0];
%  
y3=[ 3     0     0
     2     0     1
     2     1     0
     1     1     1
     1     0     2
     1     2     0
     0     2     1
     0     1     2
     0     0     3
     0     3     0];
% 
% 
y4=[ 4     0     0
     3     0     1
     3     1     0
     2     1     1
     2     0     2
     2     2     0
     1     2     1
     1     1     2
     1     0     3
     1     3     0
     0     3     1
     0     2     2
     0     1     3
     0     0     4
     0     4     0];
% 
% 
y5=[ 5     0     0 
     4     0     1
     4     1     0
     3     1     1
     3     0     2
     3     2     0
     2     2     1
     2     1     2
     2     0     3
     2     3     0
     1     3     1
     1     2     2
     1     1     3
     1     0     4
     1     4     0
     0     4     1
     0     3     2
     0     2     3
     0     1     4
     0     0     5
     0     5     0];
% 
y6=[6     0     0
     5     0     1
     5     1     0
     4     1     1
     4     0     2
     4     2     0
     3     2     1
     3     1     2
     3     0     3
     3     3     0
     2     3     1
     2     2     2
     2     1     3
     2     0     4
     2     4     0
     1     4     1
     1     3     2
     1     2     3
     1     1     4
     1     0     5
     1     5     0
     0     5     1
     0     4     2
     0     3     3
     0     2     4
     0     1     5
     0     0     6
     0     6     0];
 
muu=MD.M1(:)';
P=[MD.M2(1),MD.M2(3),MD.M2(2);MD.M2(3),MD.M2(6),MD.M2(4);MD.M2(2),MD.M2(4),MD.M2(5)];
P=P-muu'*muu;
P=P*2;
MDt=AffineTransMoms(MD,6,-sqrtm(inv(P))*muu',sqrtm(inv(P)));

% keyboard
if Nmoms==4
y=[zeros(1,length(muu));y1;y2;y3;y4];
M=[1;MDt.M1(:);MDt.M2(:);MDt.M3(:);MDt.M4(:)];
[Y,lam,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');
end

if Nmoms==5
y=[zeros(1,length(muu));y1;y2;y3;y4;y5];
M=[1;MDt.M1(:);MDt.M2(:);MDt.M3(:);MDt.M4(:);MDt.M5(:)];
[Y,lam,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');
end

if Nmoms==6
y=[zeros(1,length(muu));y1;y2;y3;y4;y5;y6([1:5],:)];
M=[1;MDt.M1(:);MDt.M2(:);MDt.M3(:);MDt.M4(:);MDt.M5(:);MDt.M6([1:5])'];
[Y,lam,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');
end

%%
if exist('CUT8PMEpdfs.mat', 'file')==2
load('CUT8PMEpdfs')
YY{end+1}=Y;
LAM{end+1}=lam;
save('CUT8PMEpdfs','YY','LAM')
else
  YY{1}=Y;
LAM{1}=lam;  
save('CUT8lamY4moms','YY','LAM')
end

dtstr=datestr(now,13);
dtstr(strfind(dtstr,':'))='-';

% keyboard
S=[MD.M2(1),MD.M2(3),MD.M2(2);MD.M2(3),MD.M2(6),MD.M2(4);MD.M2(2),MD.M2(4),MD.M2(5)];
S=S-muu'*muu;
%marg out z
[xx,yy]=meshgrid(linspace(muu(1)-3*sqrt(S(1,1)),muu(1)+3*sqrt(S(1,1)),100),linspace(muu(2)-3*sqrt(S(2,2)),muu(2)+3*sqrt(S(2,2)),100));
pent=zeros(size(xx));


[XX,W] = GH_points(muu(3),S(3,3),35);
W=W./mvnpdf(XX, muu(3),S(3,3));

iP=sqrtm(inv(P));
 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
        sp=(iP*([xx(i,j),yy(i,j),XX(k,:)]-muu)')';    
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y)*det(iP);
         end
    end
 end
figure(1)
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 hold on
 contour(xx,yy,pent,15)
 xlabel('x')
 ylabel('y')
 plot_prop_paper
 saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'fig')
  
 figure(2)
 mesh(xx,yy,pent)
 hold on
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k--','linewidth',3)
 xlabel('x')
 ylabel('y')
  plot_prop_paper
 saveas(gcf,strcat('MISSxy2_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSxy2_',num2str(kkk)),'fig')
  
  %% y-z 

[yy,zz]=meshgrid(linspace(muu(2)-3*sqrt(S(2,2)),muu(2)+3*sqrt(S(2,2)),100),linspace(muu(3)-5*sqrt(S(3,3)),muu(3)+5*sqrt(S(3,3)),100));
pent=zeros(size(yy));


[XX,W] = GH_points(muu(1),S(1,1),35);
W=W./mvnpdf(XX, muu(1),S(1,1));


 for i=1:1:size(yy,1)
    for j=1:1:size(yy,2)
        for k=1:1:length(W)
        sp=(iP*([XX(k,:),yy(i,j),zz(i,j)]-muu)')';    
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y)*det(iP);
         end
    end
 end
figure(3)
pp=sort(min(pent));
 contour(yy,zz,pent,linspace(pp(1)*1,max(max(pent)),25))
 hold on
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
  xlabel('y')
 ylabel('z')
  plot_prop_paper
  saveas(gcf,strcat('MISSyz1_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSyz1_',num2str(kkk)),'fig')
 
 figure(4)
 mesh(yy,zz,pent)
 hold on
  plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 xlabel('y')
 ylabel('z')
  plot_prop_paper
 saveas(gcf,strcat('MISSyz2_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSyz2_',num2str(kkk)),'fig')
 
 %% x-z 

[xx,zz]=meshgrid(linspace(muu(1)-3*sqrt(S(1,1)),muu(1)+3*sqrt(S(1,1)),100),linspace(muu(3)-5*sqrt(S(3,3)),muu(3)+5*sqrt(S(3,3)),100));
pent=zeros(size(xx));


[XX,W] = GH_points(muu(2),S(2,2),35);
W=W./mvnpdf(XX, muu(2),S(2,2));


 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(W)
        sp=(iP*([xx(i,j),XX(k,:),zz(i,j)]-muu)')';    
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y)*det(iP);
         end
    end
 end
figure(5)
pp=sort(min(pent));
 contour(xx,zz,pent,linspace(pp(2)*1,max(max(pent)),25))
 hold on
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 xlabel('x')
 ylabel('z')
  plot_prop_paper
 saveas(gcf,strcat('MISSxz1_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSxz1_',num2str(kkk)),'fig')
 
 figure(6)
 mesh(xx,zz,pent)
  hold on
  plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 xlabel('x')
 ylabel('z')
  plot_prop_paper
   saveas(gcf,strcat('MISSxz2_',dtstr),'pdf')
 saveas(gcf,strcat('MISSxz2_',dtstr),'fig')
 
 pause(1)
 close all
%% collision probability computation
prob_coll=zeros(1,length(R));
iP=sqrtm(inv(P));
ss=1;
for r=R
xu=r*ones(1,length(muu));
xl=-r*ones(1,length(muu));
[X,W] = GLeg_pts(35*ones(1,length(muu)),xl ,xu );
W=prod(abs(xu-xl))*W;
for k=1:1:length(W)
    XX=(iP*(X(k,:)-muu)')';
    prob_coll(ss)=prob_coll(ss)+W(k)*pdf_MaxEnt(XX,lam,Y)*det(iP);
end
ss=ss+1;
end
toc
end