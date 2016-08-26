function p=pdf_MaxEnt1D(x,lam,y)
x=x(:);

lam=lam(:);

nm=length(lam);
p=zeros(length(x),1);
for i=1:1:length(x)
p(i)=exp(sum(prod(repmat(x(i),nm,1).^y,2).*lam));
end
p=p';
end