function YY=PMEmissPDFS(MD,Nmoms,kkk,X1,X2,w1,w2,meth,filename)
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


if strcmp(meth,'MD')==1
    muu=MD.M1(:)';
    P=[MD.M2(1),MD.M2(3),MD.M2(2);MD.M2(3),MD.M2(6),MD.M2(4);MD.M2(2),MD.M2(4),MD.M2(5)];
    P=P-muu'*muu;
    P=P*2;
    MDt=AffineTransMoms(MD,Nmoms,-sqrtm(inv(P))*muu',sqrtm(inv(P)));
    
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
    
else
    % use the points  to do the miss momsnt and shifting and momoents
    n=size(X1,2);
    NN=length(w1)*length(w2);
    W=zeros(NN,1);
    X=zeros(NN,n);
    k=1;
    for i=1:1:length(w1)
        for j=1:1:length(w2)
            X(k,:) =X1(i,:)-X2(j,:);
            W(k)=w1(i)*w2(j);
            k=k+1;
        end
    end
    
    [muu,P]=MeanCov(X,W);
    % ==============  SCALING =========================
    P=P*2;
    % ==============  ------- =========================
    XX=zeros(NN,n);
    iP=sqrtm(inv(P));
    for i=1:1:NN
        XX(i,:)=iP*(X(i,:)'-muu);
    end
    
    
    if Nmoms==4
        [y1,M1]=Cal_moments_samples(XX,W,1,'raw');
        [y2,M2]=Cal_moments_samples(XX,W,2,'raw');
        [y3,M3]=Cal_moments_samples(XX,W,3,'raw');
        [y4,M4]=Cal_moments_samples(XX,W,4,'raw');
        
        y=[zeros(1,length(muu));y1;y2;y3;y4];
        M=[1;M1;M2;M3;M4];
        [Y,lam,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');
    end
    
    if Nmoms==5
        [y1,M1]=Cal_moments_samples(XX,W,1,'raw');
        [y2,M2]=Cal_moments_samples(XX,W,2,'raw');
        [y3,M3]=Cal_moments_samples(XX,W,3,'raw');
        [y4,M4]=Cal_moments_samples(XX,W,4,'raw');
        [y5,M5]=Cal_moments_samples(XX,W,4,'raw');
        
        y=[zeros(1,length(muu));y1;y2;y3;y4;y5];
        M=[1;M1;M2;M3;M4;M5];
        [Y,lam,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');
    end
    
    if Nmoms==6
        y=[zeros(1,length(muu));y1;y2;y3;y4;y5;y6([1:5],:)];
        M=[1;MDt.M1(:);MDt.M2(:);MDt.M3(:);MDt.M4(:);MDt.M5(:);MDt.M6([1:5])'];
        [Y,lam,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');
    end
    
end
% keyboard


%% collision probability computation
% prob_coll=zeros(1,length(R));
% iP=sqrtm(inv(P));
% ss=1;
% for r=R
% xu=r*ones(1,length(muu));
% xl=-r*ones(1,length(muu));
% [X,W] = GLeg_pts(35*ones(1,length(muu)),xl ,xu );
% W=prod(abs(xu-xl))*W;
% for k=1:1:length(W)
%     XX=(iP*(X(k,:)-muu)')';
%     prob_coll(ss)=prob_coll(ss)+W(k)*pdf_MaxEnt(XX,lam,Y)*det(iP);
% end
% ss=ss+1;
% end


%% Saving Data
%
% if exist(strcat(filename,'.mat'))==0
%     t=40:70;
%     YY=cell(length(t),5);
%     save(filename,'YY')
% end

iP=sqrtm(inv(P));

% pause(0.5*kkk)
%
% load(filename)
% YY{kkk,1}=Y;
% YY{kkk,2}=lam;
% YY{kkk,3}=muu;
% YY{kkk,4}=iP;
% YY{kkk,5}=@(X)pdf_MaxEnt((iP*(X(:)'-muu)')',lam,Y)*det(iP);
YY.Y=Y;
YY.lam=lam;
YY.muu=muu;
YY.iP=iP;
YY.pdf=@(X)pdf_MaxEnt((iP*(X(:)'-muu)')',lam,Y)*det(iP);
YY.Pt=P/2;
% save(filename,'YY')
% else
%     YY=cell(1,5);
% YY{1,1}=Y;
% YY{1,2}=lam;
% YY{1,3}=muu;
% YY{1,4}=iP;
% YY{1,5}=@(X)pdf_MaxEnt((iP*(X(:)'-muu)')',lam,Y)*det(iP);
% save('CUT8PMEpdfs','YY')
% end


toc

%%

return
keyboard


%%
[~,S]=MeanCov(X,W);
[xx,yy]=meshgrid(linspace(muu(1)-3*sqrt(S(1,1)),muu(1)+3*sqrt(S(1,1)),100),linspace(muu(2)-3*sqrt(S(2,2)),muu(2)+3*sqrt(S(2,2)),100));
pent=zeros(size(xx));

[Xi,Wi] = GH_points(muu(3),1*S(3,3),35);
Wi=Wi./mvnpdf(Xi, muu(3),1*S(3,3));


 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(Wi)
        sp=(iP*([xx(i,j),yy(i,j),Xi(k,:)]-muu')')';    
        pent(i,j)=pent(i,j)+Wi(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y)*det(iP);
         end
    end
 end
figure(1)
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 contour(xx,yy,pent,15)
 xlabel('x')
 ylabel('y')
%  plot_prop_paper
%  saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'fig')
  
 figure(2)
 mesh(xx,yy,pent)

 xlabel('x')
 ylabel('y')
 
%%
[~,S]=MeanCov(XX,W);


[xx,yy]=meshgrid(linspace(-3*sqrt(S(1,1)),+3*sqrt(S(1,1)),100),linspace(-3*sqrt(S(2,2)),+3*sqrt(S(2,2)),100));
pent=zeros(size(xx));

[Xi,Wi] = GH_points(0,S(3,3),35);
Wi=Wi./mvnpdf(Xi,0,S(3,3));


 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(Wi)
        sp=[xx(i,j),yy(i,j),Xi(k,:)];    
        pent(i,j)=pent(i,j)+Wi(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y);
         end
    end
 end
figure(3)
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 contour(xx,yy,pent,15)
 xlabel('x')
 ylabel('y')
%  plot_prop_paper
%  saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSxy1_',num2str(kkk)),'fig')
  
 figure(4)
 mesh(xx,yy,pent)

 xlabel('x')
 ylabel('y')
 
 %================================================


[yy,zz]=meshgrid(linspace(0-5*sqrt(S(2,2)),0+5*sqrt(S(2,2)),100),linspace(0-5*sqrt(S(3,3)),0+5*sqrt(S(3,3)),100));
pent=zeros(size(yy));

[Xi,Wi] = GH_points(0,S(1,1),35);
Wi=Wi./mvnpdf(Xi,0,S(1,1));


 for i=1:1:size(yy,1)
    for j=1:1:size(yy,2)
        for k=1:1:length(Wi)
%         sp=(iPs*([XX(k,:),yy(i,j),zz(i,j)]-muu)')';    
%         pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],YY.lam,YY.Y)*det(iP);
          sp=[Xi(k,:),yy(i,j),zz(i,j)];
        pent(i,j)=pent(i,j)+Wi(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y);
         end
    end
 end
figure(5)
pp=sort(min(pent));
 contour(yy,zz,pent,linspace(pp(1)*1,max(max(pent)),25))
%  hold on
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
  xlabel('y')
 ylabel('z')
%   plot_prop_paper
%   saveas(gcf,strcat('MISSyz1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSyz1_',num2str(kkk)),'fig')
 
 figure(6)
mesh(yy,zz,pent)
%  hold on
%   plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 xlabel('y')
 ylabel('z')
%   plot_prop_paper
%  saveas(gcf,strcat('MISSyz2_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSyz2_',num2str(kkk)),'fig')
 


[xx,zz]=meshgrid(linspace(0-3*sqrt(S(1,1)),0+3*sqrt(S(1,1)),100),linspace(0-3*sqrt(S(3,3)),0+3*sqrt(S(3,3)),100));
pent=zeros(size(xx));


[Xi,Wi] = GH_points(0,S(2,2),35);
Wi=Wi./mvnpdf(Xi, 0,S(2,2));


 for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        for k=1:1:length(Wi)
        sp=[xx(i,j),Xi(k,:),zz(i,j)];    
        pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([sp(1),sp(2),sp(3)],lam,Y);
         end
    end
 end
figure(7)
pp=sort(min(pent));
 contour(xx,zz,pent,linspace(pp(2)*1,max(max(pent)),25))
%  hold on
%  plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 xlabel('x')
 ylabel('z')
%   plot_prop_paper
%  saveas(gcf,strcat('MISSxz1_',num2str(kkk)),'pdf')
%  saveas(gcf,strcat('MISSxz1_',num2str(kkk)),'fig')
 
 figure(8)
 mesh(xx,zz,pent)
%   hold on
%   plot([R ,R ,-R ,-R ,R ],[R ,-R ,-R ,R ,R ],'k','linewidth',3)
 xlabel('x')
 ylabel('z')
%   plot_prop_paper
%    saveas(gcf,strcat('MISSxz2_',dtstr),'pdf')
%  saveas(gcf,strcat('MISSxz2_',dtstr),'fig')
%  
%  pause(1)
%  close all
end