function ceq=maxentFSOLVE(lam,y,M,X,W)

nm=size(y,1);
nq=length(W);
ceq=zeros(nm,nq);
parfor i=1:length(W)
    ceq(:,i)=W(i)*prod(repmat(X(i,:),nm,1).^y,2)*pdf_MaxEnt(X(i,:),lam,y);
    
end
ceq=sum(ceq,2)-M;


end