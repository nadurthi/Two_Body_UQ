function ProbcumMC=MC_cummCollProb(nn)

N=100000;
R=0.4;
t=40:70;
ProbcumMC=zeros(length(t),1);

SN=1:20;
SN(SN==7)=[];
for i=SN
 
    i
    rng shuffle
    reind=randperm(N);
    
    str=strcat('CollpaperSimsMC_',num2str(i));
    load(str)
    
    
    
    
    for j=1:1:N  % go throush each combo of trajectories
        
        X11mc=X1mc(t,1:3,j);
        X22mc=X2mc(t,1:3,j);
        
        D=X11mc-X22mc;
        Dnorm=sqrt(sum(D.^2,2));
        
        for tt=1:1:length(t)
            if min(Dnorm(1:tt))<=R
                ProbcumMC(tt)=ProbcumMC(tt)+1;
            end
        end
        
    end
    
end


ProbcumMC=ProbcumMC/(length(SN)*N);
save(strcat('MC_cumm_Prob_Coll_',num2str(nn)),'ProbcumMC','SN','R','t','N')