function probgh=prob_coll_norms(Xx,wx,Xy,wy,R)

d=zeros(size(Xx,1)*size(Xy,1),1);
W=zeros(size(Xx,1)*size(Xy,1),1);
kk=1;
ind=[1:1:size(Xy,1)]';
% MM1=0;
% MM2=0;
% MM3=0;
% MM4=0;
for i=1:1:size(Xx,1)
%     for j=1:1:size(Xy,1)
%         MM1=MM1+wx(i)*wy(j)*norm(Xx(i,:)-Xy(j,:))^1;
%         MM2=MM2+wx(i)*wy(j)*norm(Xx(i,:)-Xy(j,:))^2;
%         MM3=MM3+wx(i)*wy(j)*norm(Xx(i,:)-Xy(j,:))^3;
%         MM4=MM4+wx(i)*wy(j)*norm(Xx(i,:)-Xy(j,:))^4;
%         d(kk)=norm(Xx(i,:)-Xy(j,:));
%         W(kk)=wx(i)*wy(j);
%         kk=kk+1; 
         d(kk:kk+length(ind)-1)=sqrt(sum((Xx-Xy(ind,:)).^2,2));  
         W(kk:kk+length(ind)-1)=wx.*wy(ind);
        kk=kk+length(ind);
        ind=circshift(ind,1);

%     end
end
P=1*max(d);
D=d/P;
% plot(d/P,zeros(1,length(d)),'bo')
 [y1,MM1]=Cal_moments_samples(D,W,1,'raw');
 [y2,MM2]=Cal_moments_samples(D,W,2,'raw');
 [y3,MM3]=Cal_moments_samples(D,W,3,'raw');
 [y4,MM4]=Cal_moments_samples(D,W,4,'raw');
 [y5,MM5]=Cal_moments_samples(D,W,5,'raw');
  [y6,MM6]=Cal_moments_samples(D,W,6,'raw');
%    [y7,MM7]=Cal_moments_samples(D,W,7,'raw');
%      [y8,MM8]=Cal_moments_samples(D,W,8,'raw');
%     [y9,MM9]=Cal_moments_samples(D,W,9,'raw');
%       [y10,MM10]=Cal_moments_samples(D,W,10,'raw');    
M=[1;MM1;MM2;MM3;MM4;MM5;MM6];
 y=[0;y1;y2;y3;y4;y5;y6];
% y=[0;1;2;3;4];
L=1;
[Y4,lam4,xl,xu]=MaxEntPdf(y,M,0*ones(1,1),L*ones(1,1),zeros(size(y,1),1),'GLag');
[X,W] = GLeg_pts(50*ones(1), 0, R/P);
W=abs(prod(0-R/P))*W;
probgh=0;
for i=1:1:length(X)
    xx=X(i);
probgh=probgh+W(i)*pdf_MaxEnt(xx,lam4,Y4); 
end



xx=0:0.01:10;
pp=zeros(size(xx));
for i=1:1:length(xx)
    pp(i)=pdf_MaxEnt(xx(i),lam4,Y4);   
end
figure
plot(xx,pp,d(1:1e3)/P,zeros(1,length(d(1:1e3)/P)),'bo')
end