
clc
clear
close all

opt = odeset('reltol',1e-12,'abstol',1e-12);
r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
T1=[0:30:172800];

P1=blkdiag(0.01,0.01,0.01,1e-7,1e-7,1e-7);

qd_pts=@conjugate_dir_gausspts_till_8moment;

MU_cut8=cell(length(T1),5);
MU_normal=cell(length(T1),3);
for i=1:1:length(T1)
    MU_cut8{i,5}=0;
    MU_normal{i,3}=0;
end

R=diag([(0.1)^2]);    
hn=1;
hx=@(x)[sqrt(x(1)^2+x(2)^2+x(3)^2)];


[t,ytruth]=ode45(@twoBody,T1,r1,opt);

rr1=mvnrnd(r1,3*P1);

[X,w]=qd_pts(rr1,P1);
[y1,M1]=Cal_moments_samples(X,w,1,'central');
[y2,M2]=Cal_moments_samples(X,w,2,'central');
[y3,M3]=Cal_moments_samples(X,w,3,'central');
[y4,M4]=Cal_moments_samples(X,w,4,'central');

MU_cut8{1,1}=M1;
MU_cut8{1,2}=M2;
MU_cut8{1,3}=M3;
MU_cut8{1,4}=M4;

     mu_normal=rr1(:);
     P_normal=P1;
     MU_normal{1,1}=mu_normal;
     MU_normal{1,2}=P_normal;
    
     Xprior=X;
     Wprior=w;
     Npts=length(Wprior);
     
     pk=1;
HIGHrmse=zeros(length(T1),1);
LOWrmse=zeros(length(T1),1);

KK=300:300:length(T1);
  for k=1:1:length(KK)
disp([num2str(KK(k)),' of ',length(T1)] )
        
        if rem(KK(k),600)==0
            
            %% common stuff
            ym=hx(ytruth(KK(k),:))+sqrtm(R)*rand(hn,1);
            
            %% Higher Moment
            parfor i=1:1:Npts
                [t,x]=ode45(@twoBody,[T1(pk),T1(KK(k))],Xprior(i,:),opt);
                Xprior(i,:)=x(end,:);
            end
            %forecast mean and cov at time k
            [mu,P]=MeanCovPts(Wprior,Xprior);
            % observation mean dn cov at time k
            Zprior=zeros(Npts,hn);
            for i=1:1:Npts
                Zprior(i,:)=hx(Xprior(i,:));
            end
            [muz,Pz]=MeanCovPts(Wprior,Zprior);
            Pz=Pz+R;
            % cross covariance 
            Pcc=0;
            for i=1:1:length(w)
                Pcc=Pcc+Wprior(i)*(Xprior(i,:)'-mu(:))*(Zprior(i,:)-muz(:)');
            end
            %kalman gain
            K=Pcc/Pz;
            %update
            mu_post=mu(:)+K*(ym-muz(:));
            P_post=P-K*Pz*K';
            
            [Xpost,Wpost]=qd_pts(mu_post,P_post);
            
             [Xpost,Wpost]=Reweightpts(4,Xprior,Wprior,hx,R,ym,Xpost,Wpost);
             Wprior=Wpost;
             Xprior=Xpost;
             [y1,M1]=Cal_moments_samples(Xprior,Wprior,1,'central');
             [y2,M2]=Cal_moments_samples(Xprior,Wprior,2,'central');
             [y3,M3]=Cal_moments_samples(Xprior,Wprior,3,'central');
             [y4,M4]=Cal_moments_samples(Xprior,Wprior,4,'central');
             
             MU_cut8{KK(k),1}=M1;
             MU_cut8{KK(k),2}=M2;
             MU_cut8{KK(k),3}=M3;
             MU_cut8{KK(k),4}=M4;

          %% Only 2 monents   
          [Xnormal,Wnormal]=qd_pts(mu_normal,P_normal);
          parfor i=1:1:length(Wnormal)
              [t,x]=ode45(@twoBody,[T1(pk),T1(KK(k))],Xnormal(i,:),opt);
              Xnormal(i,:)=x(end,:);
          end
          [mu_normal,P_normal]=MeanCovPts(Wnormal,Xnormal);
          Znormal=zeros(length(Wnormal),hn);
          for i=1:1:length(Wnormal)
              Znormal(i,:)=hx(Xnormal(i,:));
          end
          [muznormal,Pznormal]=MeanCovPts(Wnormal,Znormal);
          Pznormal=Pznormal+R;
          % cross covariance
          Pccnormal=0;
          for i=1:1:length(w)
              Pccnormal=Pccnormal+Wnormal(i)*(Xnormal(i,:)'-mu_normal(:))*(Znormal(i,:)-muznormal(:)');
          end
          %kalman gain
          K=Pccnormal/Pznormal;
          %update
          mu_normal=mu_normal(:)+K*(ym-muznormal(:));
          P_normal=P_normal-K*Pznormal*K';
          
          MU_normal{KK(k),1}=mu_normal;
          MU_normal{KK(k),2}=P_normal;
          
          %% changing prev update time step
          pk=KK(k);
          
          
        else
             %% Higher monents  
%              keyboard
             XX=zeros(size(Xprior));
            parfor i=1:1:Npts
                [t,x]=ode45(@twoBody,[T1(pk),T1(KK(k))],Xprior(i,:),opt);
                XX(i,:)=x(end,:);
            end
            
            [y1,M1]=Cal_moments_samples(XX,Wprior,1,'raw');
            [y2,M2]=Cal_moments_samples(XX,Wprior,2,'raw');
            [y3,M3]=Cal_moments_samples(XX,Wprior,3,'raw');
            [y4,M4]=Cal_moments_samples(XX,Wprior,4,'raw');
            MU_cut8{KK(k),1}=M1;
            MU_cut8{KK(k),2}=M2;
            MU_cut8{KK(k),3}=M3;
            MU_cut8{KK(k),4}=M4;
             %% Only 2 monents  
             [Xnormal,Wnormal]=qd_pts(mu_normal,P_normal);
             XX=zeros(size(Xnormal));
             parfor i=1:1:length(Wnormal)
                 [t,x]=ode45(@twoBody,[T1(pk),T1(KK(k))],Xnormal(i,:),opt);
                 XX(i,:)=x(end,:);
             end
             [mu_inter,P_inter]=MeanCovPts(Wnormal,XX);
             MU_normal{KK(k),1}=mu_inter;
             MU_normal{KK(k),2}=P_inter;
        end

        muHIGH=zeros(length(KK(1:k)),6);
        muLOW= zeros(length(KK(1:k)),6);
        KKK=KK(1:k);
        for ii=1:1:length(KK(1:k))
           muHIGH(ii,:)=  MU_cut8{KKK(ii),1}(1:6);
           muLOW(ii,:)= MU_normal{KKK(ii),1}(1:6);
        end
        
        HIGHrmse(KK(k))=RMSEMat(ytruth(KK(1:k),:),muHIGH);
        LOWrmse(KK(k))=RMSEMat(ytruth(KK(1:k),:),muLOW);
        
        plot(T1(KK(1:k))/(60*60),HIGHrmse(KK(1:k)),'ks-',T1(KK(1:k))/(60*60),LOWrmse(KK(1:k)),'b*-','MarkerSize',6)
        legend('HigherMoms','LowMoms')
%       plot3(ytruth(1:k,1),ytruth(1:k,2),ytruth(1:k,3),'ro-',muHIGH(1:k,1),muHIGH(1:k,2),muHIGH(1:k,3),'ks-',muLOW(1:k,1),muLOW(1:k,2),muLOW(1:k,3),'b*-','MarkerSize',6,'linewidth',2) 
%         legend('truth','HigherMoms','LowMoms')
        pause(0.5)
%         keyboard
  end
 

  


