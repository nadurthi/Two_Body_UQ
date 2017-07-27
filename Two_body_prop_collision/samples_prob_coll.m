function prob_collmc=samples_prob_coll(X1,X2,R,method)
dim=size(X1,2);
N=size(X1,1);
D=zeros(N,dim);
ind=randperm(N);
prob_collmc=0;
for i=1:1:N
    D(i,:)=(X1(ind(i),:)-X2(i,:));
    if method==1
    if sum(sign(R*ones(1,dim)-D(i,:)))==dim && sum(sign(R*ones(1,dim)+D(i,:)))==dim
    prob_collmc=prob_collmc+1;
    end
    elseif method==2
    if norm(D(i,:))<R
    prob_collmc=prob_collmc+1;
    end
    end
    
end
prob_collmc=prob_collmc/size(D,1);