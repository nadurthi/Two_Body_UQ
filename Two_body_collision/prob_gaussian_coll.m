function prob=prob_gaussian_coll(mu1,P1,mu2,P2,R,kkk)
% [Xx,wx]=GH_points(mu1(:),P1,10);
% [Xy,wy]=GH_points(mu2(:),P2,10);
% 
% prob_full=prob_coll_nongauss(Xx,wx,Xy,wy,R);
% probnorm=prob_coll_norms(Xx,wx,Xy,wy,R);

 muc=mu1-mu2;
 Pc=P1+P2;
% [Xd,wd]=GH_points(muc(:),Pc,10);
%  [y1,M1d]=Cal_moments_samples(Xd ,wd,1,'raw');
%  [y2,M2d]=Cal_moments_samples(Xd ,wd,2,'raw');
%  [y3,M3d]=Cal_moments_samples(Xd ,wd,3,'raw');
%  [y4,M4d]=Cal_moments_samples(Xd ,wd,4,'raw');
% keyboard
 dim=length(muc);
 
 %% x-y
 pp=[1,2];
 [xx,yy]=meshgrid(linspace(muc(pp(1))-3*sqrt(Pc(pp(1),pp(1))),muc(pp(1))+3*sqrt(Pc(pp(1),pp(1))),100),linspace(muc(pp(2))-3*sqrt(Pc(pp(2),pp(2))),muc(pp(2))+3*sqrt(Pc(pp(2),pp(2))),100));
 pent=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
 pent(i,j)= mvnpdf([xx(i,j),yy(i,j)],muc(pp),Pc(pp,pp));
    end
end

figure(1)
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 hold on
 contour(xx,yy,pent,15)
 xlabel('x')
 ylabel('y')
 plot_prop_paper
 saveas(gcf,strcat('MISSGaussxy1_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSGaussxy1_',num2str(kkk)),'fig')
  
 figure(2)
 mesh(xx,yy,pent)
 view(153,58)
 hold on
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k--','linewidth',3)
 xlabel('x')
 ylabel('y')
  plot_prop_paper
 saveas(gcf,strcat('MISSGaussxy2_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSGaussxy2_',num2str(kkk)),'fig')
 
 %% y-z
 pp=[2,3];
 [xx,yy]=meshgrid(linspace(muc(pp(1))-3*sqrt(Pc(pp(1),pp(1))),muc(pp(1))+3*sqrt(Pc(pp(1),pp(1))),100),linspace(muc(pp(2))-3*sqrt(Pc(pp(2),pp(2))),muc(pp(2))+3*sqrt(Pc(pp(2),pp(2))),100));
 pent=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
 pent(i,j)= mvnpdf([xx(i,j),yy(i,j)],muc(pp),Pc(pp,pp));
    end
end

figure(3)
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 hold on
 contour(xx,yy,pent,15)
 xlabel('y')
 ylabel('z')
 plot_prop_paper
 saveas(gcf,strcat('MISSGaussyz1_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSGaussyz1_',num2str(kkk)),'fig')
  
 figure(4)
 mesh(xx,yy,pent)
 view(153,58)
 hold on
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k--','linewidth',3)
 xlabel('y')
 ylabel('z')
  plot_prop_paper
 saveas(gcf,strcat('MISSGaussyz2_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSGaussyz2_',num2str(kkk)),'fig')
  %% x-z
 pp=[1,3];
 [xx,yy]=meshgrid(linspace(muc(pp(1))-3*sqrt(Pc(pp(1),pp(1))),muc(pp(1))+3*sqrt(Pc(pp(1),pp(1))),100),linspace(muc(pp(2))-3*sqrt(Pc(pp(2),pp(2))),muc(pp(2))+3*sqrt(Pc(pp(2),pp(2))),100));
 pent=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
 pent(i,j)= mvnpdf([xx(i,j),yy(i,j)],muc(pp),Pc(pp,pp));
    end
end

figure(5)
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k','linewidth',3)
 hold on
 contour(xx,yy,pent,15)
 xlabel('x')
 ylabel('z')
 plot_prop_paper
 saveas(gcf,strcat('MISSGaussxz1_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSGaussxz1_',num2str(kkk)),'fig')
  
 figure(6)
 mesh(xx,yy,pent)
 view(153,58)
 hold on
 plot([R(13),R(13),-R(13),-R(13),R(13)],[R(13),-R(13),-R(13),R(13),R(13)],'k--','linewidth',3)
 xlabel('x')
 ylabel('z')
  plot_prop_paper
 saveas(gcf,strcat('MISSGaussxz2_',num2str(kkk)),'pdf')
 saveas(gcf,strcat('MISSGaussxz2_',num2str(kkk)),'fig')
 
 close all
%%
 probanal=zeros(1,length(R));
 
 prob=probanal;
 return
 
 ss=1;
 for r=R
 xl=-r*ones(1,dim);
 xu=r*ones(1,dim);
 [xint,wint] = GLeg_pts(35*ones(1,dim),xl, xu);
 wint=prod(abs(xu-xl))*wint;
  for i=1:1:length(wint)
    probanal(ss)=probanal(ss)+wint(i)*mvnpdf(xint(i,:),muc,Pc); 
  end
 ss=ss+1;
 end

  
% 
% prob_mc=zeros(1,10);
% for j=1:1:10
% Xmc1=mvnrnd(mu1,P1,1e6);
% Xmc2=mvnrnd(mu2,P2,1e6);
% N=size(Xmc1,1);
% D=zeros(N,dim);
% ind=randperm(N);
% prob_mc(j)=0;
% for i=1:1:N
%     D(i,:)=(Xmc1(ind(i),:)-Xmc2(i,:));
%     if sum(sign(R*ones(1,dim)-D(i,:)))==dim && sum(sign(R*ones(1,dim)+D(i,:)))==dim
%     prob_mc(j)=prob_mc(j)+1;
%     end
% end
% j
% prob_mc(j)=prob_mc(j)/N;
% end
% figure
% plot3(Xmc1(:,1),Xmc1(:,2),Xmc1(:,3),'bo',Xmc2(:,1),Xmc2(:,2),Xmc2(:,3),'ro')
prob=probanal;
%   prob=[prob_full,probnorm,probanal,max(prob_mc)];
end