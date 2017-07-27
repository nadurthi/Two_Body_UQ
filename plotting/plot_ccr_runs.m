k=721
XsigSat=XsigSat_ut;
nnn=length(XsigSat{1,1}(1,:));
nsat=size(XsigSat,1);
% nsat=10;
H=zeros(k,nsat);
s={};
err=zeros(k,nsat);
NN=[1:10,12:21,23:25]
nnn=nnn-2;
% NN=1:25
for i=NN
for t=1:1:k
err(t,i)=sqrt(1/t*sum(sqrt(1/nnn*sum((XsigSat{i,1}(1:t,:)-ytruth{i}(1:t,:)).^2)).^2));
end
end
HHut=err;
XsigSat=XsigSat_cut8;
% nnn=length(XsigSat{1,1}(1,:));
% nsat=size(XsigSat,1);
H=zeros(k,nsat);
s={};
err=zeros(k,nsat);
for i=NN
for t=1:1:k
err(t,i)=sqrt(1/t*sum(sqrt(1/nnn*sum((XsigSat{i,1}(1:t,:)-ytruth{i}(1:t,:)).^2)).^2));
end
end
HHcut8=err;

for pp=1:1:25
figure
plot(Tvec/3600,HHcut8(:,pp),Tvec/3600,HHut(:,pp),'--','linewidth',2)
legend('CUT8','UT')
end

plot(Tvec,Hcut8err,Tvec,Huterr,'--')  
legend('cut8','ut')      


pp=1
eigcut=0;
eigu=0;
for ss=1:1:k
    eigcut(ss)=max(eig(reshape(XsigSat_cut8{pp,2}(ss,:),6,6)));
    eigu(ss)=max(eig(reshape(XsigSat_ut{pp,2}(ss,:),6,6)));
end   
figure   
 plot(Tvec/3600,eigcut,Tvec/3600,eigu,'--')   
 
 figure
 plot(Tvec/3600,XsigSat_cut8{pp,1}(:,4),Tvec/3600,XsigSat_ut{pp,1}(:,4),'--')    
 
 
 Radmodel.R=@(Srad)reshape(MeasCov(Srad,:),Radmodel.hn,Radmodel.hn);
Radmodel.R(1)     

flag=[];
for k=1:1:289
try
    S=MeasPairs_cut8{k}-MeasPairs_ut{k};
catch
   flag=vertcat(flag,k); 
   k
end
if sum(sum(abs(S)))>0
    flag=vertcat(flag,k);
end
end 

%%
MSup=min_time_btw_updates(nsat,289,MeasPairs_ut)
MM=zeros(20,111);
for i=1:1:nsat
    for j=1:1:length(MSup{i})
MM(i,j)=MSup{i}(j);
    end        
end     
xlswrite('MeasTintervals',MM,'A1')
                               


mu=398601.2;
[a,e,E,w,i,Om] = XYZ2OE(Xsat0(1,1:3)',Xsat0(1,4:6)',mu) 

TP=pi*sqrt(a^3/mu)





for ct=1:2000 % Index for number of RSOs considered.
    P(1,1,ct)=Covariance(6*ct-5,1);
    P(2,2,ct)=Covariance(6*ct-4,1);
    P(3,3,ct)=Covariance(6*ct-3,1);
    P(1,2,ct)=Covariance(6*ct-2,1);
    P(1,3,ct)=Covariance(6*ct-1,1);
    P(2,3,ct)=Covariance(6*ct,1);
    P(2,1,ct)=P(1,2,ct);
    P(3,1,ct)=P(1,3,ct);
    P(3,2,ct)=P(2,3,ct):
end

Xsat1=zeros(2000,6);

disp(1)          
for i=1:1:2000  
    Xsat1(i,1:3)=pos_store(1,(i-1)*3+1:i*3);
    Xsat1(i,4:6)=vel_store(1,(i-1)*3+1:i*3);
end
disp(2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       