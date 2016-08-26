function [mu,P]=ptswts2muP(X,w)
w=w(:);
dim=size(X,2);
N=size(X,1);
mu=sum(repmat(w,1,dim).*X,1);
X=X-repmat(mu,N,1);
P=zeros(dim,dim);
for i=1:1:N
    P=P+w(i)*X(i,:)'*X(i,:);
end

end
