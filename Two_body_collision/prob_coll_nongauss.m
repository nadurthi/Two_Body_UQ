function prob_coll=prob_coll_nongauss(Xx,wx,Xy,wy,R)
dim=size(Xx,2);
%% first calculate the moments of the difference
%  [y1,M1x]=Cal_moments_samples(Xx ,wx,1,'raw');
%  [y2,M2x]=Cal_moments_samples(Xx ,wx,2,'raw');
%  [y3,M3x]=Cal_moments_samples(Xx ,wx,3,'raw');
%  [y4,M4x]=Cal_moments_samples(Xx ,wx,4,'raw');
%  Mx={M1x';M2x';M3x';M4x'};
%  
%  [y1,M1y]=Cal_moments_samples(Xy ,wy,1,'raw');
%  [y2,M2y]=Cal_moments_samples(Xy ,wy,2,'raw');
%  [y3,M3y]=Cal_moments_samples(Xy ,wy,3,'raw');
%  [y4,M4y]=Cal_moments_samples(Xy ,wy,4,'raw');
%  My={M1y';M2y';M3y';M4y'};
% 
%  MD=Diff_rv_moments(Mx,My);

%  tic
 D=zeros(size(Xx,1)*size(Xx,1),3);
 W=zeros(size(D,1),1);
 k=1;
 for i=1:1:size(Xx,1)
     for j=1:1:size(Xy,1)
     D(k,:)=Xx(i,:)-Xy(j,:);
     W(k)=wx(i)*wy(j);
     k=k+1;
     end
 end

% %  MoMs=missdistpts(Xx,wx,Xy,wy);
%   [y1,M1x]=Cal_moments_samples(D ,W,1,'raw');
%  [y2,M2x]=Cal_moments_samples(D ,W,2,'raw');
%  [y3,M3x]=Cal_moments_samples(D ,W,3,'raw');
%  [y4,M4x]=Cal_moments_samples(D ,W,4,'raw');
%   [y5,M5x]=Cal_moments_samples(D ,W,5,'raw');
%    [y6,M6x]=Cal_moments_samples(D ,W,6,'raw');
%  MMx.M1=M1x;
%  MMx.M2=M2x;
%  MMx.M3=M3x;
%  MMx.M4=M4x;
%  MMx.M5=M5x;
%  MMx.M6=M6x;
%  
%  MDth={M1Dth';M2Dth';M3Dth';M4Dth'};
% keyboard
 [muu2,P2]=ptswts2muP(D,W);
Dt=zeros(size(D));
for i=1:1:size(D,1)
   Dt(i,:)=sqrtm(inv(P2))*(D(i,:)-muu2)'; 
end
[y1,M1tDth]=Cal_moments_samples(Dt ,W,1,'raw');
 [y2,M2tDth]=Cal_moments_samples(Dt ,W,2,'raw');
 [y3,M3tDth]=Cal_moments_samples(Dt ,W,3,'raw');
 [y4,M4tDth]=Cal_moments_samples(Dt ,W,4,'raw');
  [y5,M5tDth]=Cal_moments_samples(Dt ,W,5,'raw');
[y6,M6tDth]=Cal_moments_samples(Dt ,W,6,'raw');
 MDt2={M1tDth';M2tDth';M3tDth';M4tDth';M5tDth';M6tDth'};
% keyboard

%% Anal methods
% [y1,M1x]=Cal_moments_samples(Xx ,wx,1,'raw');
% [y2,M2x]=Cal_moments_samples(Xx ,wx,2,'raw');
% [y3,M3x]=Cal_moments_samples(Xx ,wx,3,'raw');
% [y4,M4x]=Cal_moments_samples(Xx ,wx,4,'raw');
% [y5,M5x]=Cal_moments_samples(Xx ,wx,5,'raw');
% [y6,M6x]=Cal_moments_samples(Xx ,wx,6,'raw');
% 
% [y1,M1y]=Cal_moments_samples(Xy ,wy,1,'raw');
% [y2,M2y]=Cal_moments_samples(Xy ,wy,2,'raw');
% [y3,M3y]=Cal_moments_samples(Xy ,wy,3,'raw');
% [y4,M4y]=Cal_moments_samples(Xy ,wy,4,'raw');
% [y5,M5y]=Cal_moments_samples(Xy ,wy,5,'raw');
% [y6,M6y]=Cal_moments_samples(Xy ,wy,6,'raw');
% 
% Mx=struct('M1',M1x,'M2',M2x,'M3',M3x,'M4',M4x,'M5',M5x,'M6',M6x);
% My=struct('M1',M1y,'M2',M2y,'M3',M3y,'M4',M4y,'M5',M5y,'M6',M6y);
% 
% MD=MissMomsFromIndividMoms(Mx,My,6);
muu=MD.M1(:)';
P=[MD.M2(1),MD.M2(3),MD.M2(2);MD.M2(3),MD.M2(6),MD.M2(4);MD.M2(2),MD.M2(4),MD.M2(5)];
P=P-muu'*muu;
MDt=AffineTransMoms(MD,6,-sqrtm(inv(P))*muu',sqrtm(inv(P)));
keyboard
%% now convert the moments to 0 mean and I cov
% muu=MD{1};
% MDs=MOMS_shift(MD,muu);
% if dim==2
% P=1*[MDs{2}(1),MDs{2}(2);MDs{2}(2),MDs{2}(3)];
% elseif dim==3    
% P=1*[MDs{2}(1),MDs{2}(3),MDs{2}(2);MDs{2}(3),MDs{2}(6),MDs{2}(4);MDs{2}(2),MDs{2}(4),MDs{2}(5)];    
% end
% MDt=MOMS_tranform(MDs,inv(sqrtm(P)));

% [MDth{1}',MD{1}']
% [MDth{2}',MD{2}']
% [MDth{3}',MD{3}']
% [MDth{4}',MD{4}']
% 
% [MDtth{1}',MDt{1}']
% [MDtth{2}',MDt{2}']
% [MDtth{3}',MDt{3}']
% [MDtth{4}',MDt{4}']
%% further scaling


%%%%%%%%%%%%%%% 2nd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y=[zeros(1,length(muu));y1;y2];
% M=[1;MDt{1}';MDt{2}'];
% [Y2,lam2,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),zeros(size(y),1));

%%%%%%%%%%%%%%% 3rd MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y=[zeros(1,length(muu));y1;y2;y3];
% lam0=[lam2;zeros(size(y,1)-length(lam2),1)];
% M=[1;MDt{1}';MDt{2}';MDt{3}'];
% [Y3,lam3,xl,xu]=MaxEntPdf(y,M,-1.5*s*ones(1,ns),1.5*s*ones(1,ns),lam0);

%%%%%%%%%%%%%%% 4th MOMS constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=[zeros(1,length(muu));y1;y2;y3;y4];
M=[1;MDt.M1(:);MDt.M2(:);MDt.M3(:);MDt.M4(:)];
% M([33:34])=[];
% y([33:34],:)=[];
[y,M]

[Y4,lam4,xl,xu]=MaxEntPdf(y,M,-1.5*ones(1,dim),1.5*ones(1,dim),zeros(size(y,1),1),'GH');


% N=size(Xmc1,1);
% D=zeros(N,dim);
% ind=randperm(N);
% prob_collmc=0;
% for i=1:1:N
%     D(i,:)=(Xmc1(ind(i),:)-Xmc2(i,:));
%     if sum(sign(R*ones(1,dim)-D(i,:)))==dim && sum(sign(R*ones(1,dim)+D(i,:)))==dim
%     prob_collmc=prob_collmc+1;
%     end
% end
% prob_collmc=prob_collmc/size(D,1);
% muD=mean(D,1);
% PD=cov(D);
% for i=1:1:N
%     D(i,:)=sqrtm(inv(PD))*(D(i,:)-muD)';
% end
% 
% [xx,yy]=meshgrid(linspace(-5,5,100),linspace(-5,5,100));
% [X,W] = GH_points(zeros(1,1),0.5*eye(1),35);
% W=W./mvnpdf(X,zeros(1,1),0.5*eye(1));
% pent=zeros(size(xx));
% 
% XX=zeros(size(xx));
% YY=zeros(size(xx));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         pent(i,j)=pdf_MaxEnt([xx(i,j),yy(i,j)],lam4,Y4)*det(sqrtm(inv(P)));
%         XXt=sqrtm(P)*[xx(i,j),yy(i,j)]'+muu';
%         XX(i,j)=XXt(1);
%         YY(i,j)=XXt(2);
% %         for k=1:1:length(W)
% %  pent(i,j)=pent(i,j)+W(k)*pdf_MaxEnt([xx(i,j),yy(i,j),X(k,:)],lam4,Y4);
% %         end
%     end
% end 
% 
%  figure
%  plot(D(:,1),D(:,2),'ro')
%   hold on
%  contour(XX,YY,pent,15,'linewidth',2)
%    xlabel('x')
%    ylabel('y')

% linspace(0.001,max(max(pent)/3),15)

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
    prob_coll(ss)=prob_coll(ss)+W(k)*pdf_MaxEnt(XX,lam4,Y4)*det(iP);
end
ss=ss+1;
end
toc
end