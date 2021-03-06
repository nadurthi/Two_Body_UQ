function MC_coll_ccr(numWorkers,ss)
Ncollmc=500;
%   matlabpool('local', numWorkers)
dim=3;
R=[0.01:0.01:0.1,0.2:0.1:0.5];
Tt=50:65;
prob_collmc_square=cell(1,Ncollmc);
prob_collmc_circle=cell(1,Ncollmc);
for varr=1:1:Ncollmc
    prob_collmc_square{varr}=zeros(length(R),length(Tt));
    prob_collmc_circle{varr}=zeros(length(R),length(Tt));
end

Nsim=1:51;
for f=Nsim
    load(strcat('eg1CollisionBodies',num2str(f)))
N=size(X1mc,3);
t=size(X1mc,1);
X1=zeros(N,3);
X2=zeros(N,3);


kt=1;
for t=Tt 
    
    for j=1:1:N
            X1(j,:)=X1mc(t,1:3,j);
            X2(j,:)=X2mc(t,1:3,j);
    end

kr=1;
for r=R
    
    for varr=1:1:Ncollmc

        ind=randperm(N);
        XX=X1(ind,:);
        D=XX-X2;
        S=sqrt(sum(D.^2,2));
        ncircle=length(find(S<=r));
        D=sum(sign(r*ones(N,dim)-D),2)+sum(sign(r*ones(N,dim)+D),2);
        nsquare=length(find(D==6));
        prob_collmc_square{varr}(kr,kt)=prob_collmc_square{varr}(kr,kt)+nsquare;
        prob_collmc_circle{varr}(kr,kt)=prob_collmc_circle{varr}(kr,kt)+ncircle;
    end %multiple varr loop
% prob_collmc(kr,kt)=prob_collmc(kr,kt)+samples_prob_coll(X1,X2,r,1);
kr=kr+1;
[f,kt,kr]
end
kt=kt+1;
end %time loop


end% load loop

f=length(Nsim)*1e5;
for varr=1:1:Ncollmc
    prob_collmc_square{varr}=prob_collmc_square{varr}/f;
    prob_collmc_circle{varr}=prob_collmc_circle{varr}/f;
end
save(strcat('MCconjProbComps_',num2str(ss)),'prob_collmc_circle','prob_collmc_square','Tt','R','f')


  % all done, close pool of workers
%   matlabpool close
exit;