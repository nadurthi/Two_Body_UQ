%% gaussian moment collision in 2d
dim=3;
mu1=[2,2,2];
mu2=[-2,2,2];

% P1=[1,-0.5;-0.5,1];
% P2=[1,0.5;0.5,1];
P1=eye(3);
P2=2*eye(3);

N=100000;
X1 = mvnrnd(mu1,P1,N); 
X2 = mvnrnd(mu2,P2,N); 

% plot(X1(:,1),X1(:,2),'bo',X2(:,1),X2(:,2),'ro')

% NN=1000;
% nn=zeros(1,NN);
% for j=1:1:NN
d=zeros(N,1);
D=zeros(N,dim);
ind=randperm(N);
% k=1;
n=0;
for i=1:1:N
%     for j=1:1:N
    d(i)=norm(X1(ind(i),:)-X2(i,:));
    D(i,:)=(X1(ind(i),:)-X2(i,:));
    if d(i)<0.5
        n=n+1;
    end
%     k=k+1;
%     end
end
% nn(j)=n;
% end
% hist(nn);

w=ones(1,N)/N;
 [y1,M1D]=Cal_moments_samples(D,w,1,'raw');
 [y2,M2D]=Cal_moments_samples(D,w,2,'raw');
 [y3,M3D]=Cal_moments_samples(D,w,3,'raw');
 [y4,M4D]=Cal_moments_samples(D,w,4,'raw');
 
 
 %% analytical difference bn Gaussian pdf is gaussian
 muc=mu1-mu2;
 Pc=P1+P2;
%  Xc = mvnrnd(muc,Pc,N); 
[Xc,wc]=GH_points(muc',Pc,5);

  [y1,M1c]=Cal_moments_samples(Xc ,wc,1,'raw');
 [y2,M2c]=Cal_moments_samples(Xc ,wc,2,'raw');
 [y3,M3c]=Cal_moments_samples(Xc ,wc,3,'raw');
 [y4,M4c]=Cal_moments_samples(Xc ,wc,4,'raw');
 
 
 %% Explicit moment cal from difference 
[X1,w1]=GH_points(mu1',P1,5);
[X2,w2]=GH_points(mu2',P2,5);
  [y1,M1]=Cal_moments_samples(X1 ,w1,1,'raw');
 [y2,M2]=Cal_moments_samples(X1 ,w1,2,'raw');
 [y3,M3]=Cal_moments_samples(X1 ,w1,3,'raw');
 [y4,M4]=Cal_moments_samples(X1 ,w1,4,'raw');
 Mx={M1';M2';M3';M4'};
 
 [y1,M1]=Cal_moments_samples(X2 ,w2,1,'raw');
 [y2,M2]=Cal_moments_samples(X2 ,w2,2,'raw');
 [y3,M3]=Cal_moments_samples(X2 ,w2,3,'raw');
 [y4,M4]=Cal_moments_samples(X2 ,w2,4,'raw');
 My={M1';M2';M3';M4'};
 
 
 MDD=Diff_rv_moments(Mx,My)
 
 
 
 %% Generating a special case 
 X1=Data_TBP_at699mc(:,1:2);
 X1cut=Data_TBP_at699cut8(:,1:2);
 wcut=wcut8;
 Ncut=size(X1cut,1);
 N=size(X1,1);
 X2=X1+repmat([5,0],N,1);
X2cut=X1cut+repmat([5,0],Ncut,1);
 mu=mean(X2,1);
 mucut=mean(X2cut,1);
 X2=X2-repmat(mu,N,1);
 X2cut=X2cut-repmat(mucut,Ncut,1);
R=[cosd(60),sind(60);-sind(60),cosd(60)];
 for i=1:1:N
     X2(i,:)=R*X2(i,:)'+mu';
     if i<=Ncut
     X2cut(i,:)=R*X2cut(i,:)'+mucut';
     end
 end
 w=ones(N,1)/N;
   [y1,M1]=Cal_moments_samples(X1 ,w,1,'raw');
 [y2,M2]=Cal_moments_samples(X1 ,w,2,'raw');
 [y3,M3]=Cal_moments_samples(X1 ,w,3,'raw');
 [y4,M4]=Cal_moments_samples(X1 ,w,4,'raw');
 Mx={M1';M2';M3';M4'};
 
 [y1,M1]=Cal_moments_samples(X2 ,w,1,'raw');
 [y2,M2]=Cal_moments_samples(X2 ,w,2,'raw');
 [y3,M3]=Cal_moments_samples(X2 ,w,3,'raw');
 [y4,M4]=Cal_moments_samples(X2 ,w,4,'raw');
 My={M1';M2';M3';M4'};
 
 [y1,M1]=Cal_moments_samples(X1cut ,wcut,1,'raw');
 [y2,M2]=Cal_moments_samples(X1cut ,wcut,2,'raw');
 [y3,M3]=Cal_moments_samples(X1cut ,wcut,3,'raw');
 [y4,M4]=Cal_moments_samples(X1cut ,wcut,4,'raw');
 Mxcut={M1';M2';M3';M4'};
 
 [y1,M1]=Cal_moments_samples(X2cut ,wcut,1,'raw');
 [y2,M2]=Cal_moments_samples(X2cut ,wcut,2,'raw');
 [y3,M3]=Cal_moments_samples(X2cut ,wcut,3,'raw');
 [y4,M4]=Cal_moments_samples(X2cut ,wcut,4,'raw');
 Mycut={M1';M2';M3';M4'};
 
 MDD=Diff_rv_moments(Mx,My);
 MDDcut=Diff_rv_moments(Mxcut,Mycut);
  
 dim=2;
 d=zeros(N,1);
D=zeros(N,dim);
ind=randperm(N);
% k=1;
n=0;
for i=1:1:N
%     for j=1:1:N
    d(i)=norm(X1(ind(i),:)-X2(i,:));
    D(i,:)=(X1(ind(i),:)-X2(i,:));
    if d(i)<0.5
        n=n+1;
    end
%     k=k+1;
%     end
end
 [y1,M1D]=Cal_moments_samples(D,w,1,'raw');
 [y2,M2D]=Cal_moments_samples(D,w,2,'raw');
 [y3,M3D]=Cal_moments_samples(D,w,3,'raw');
 [y4,M4D]=Cal_moments_samples(D,w,4,'raw');
 
 [M1D,MDD{1}',MDDcut{1}']
 [M2D,MDD{2}',MDDcut{2}']
 [M3D,MDD{3}',MDDcut{3}']
 [M4D,MDD{4}',MDDcut{4}']
 
 
 %% test shift cal moms
 A=eye(3);
 Mx={M1D';M2D';M3D';M4D'};
 MDD=MOMS_tranform(Mx,A);
 
 [Mx{1}',MDD{1}']
 [Mx{2}',MDD{2}']
 [Mx{3}',MDD{3}']
 [Mx{4}',MDD{4}']