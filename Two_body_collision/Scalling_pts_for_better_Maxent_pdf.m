Xmc=Data_TBP_at699mc(:,1:2);
plot(Xmc(:,1),Xmc(:,2),'bo')
P=cov(Xmc);
mu=mean(Xmc,1)';
sqP=sqrtm(inv(P));

Ymc=zeros(size(Xmc));
for i=1:1:length(Xmc)
Ymc(i,:)=sqP*(Xmc(i,:)'-mu);
end
figure
plot(Ymc(:,1),Ymc(:,2),'bo')

bl=min(Xmc,[],1)-(1/5)*abs(min(Xmc,[],1)-mu');
bu=max(Xmc,[],1)+(1/5)*abs(max(Xmc,[],1)-mu');
[Zmc,d]=transform_domain(Xmc,bl,bu,-ones(1,size(Xmc,2)),ones(1,size(Xmc,2)));
figure
plot(Zmc(:,1),Zmc(:,2),'bo')
axis([-1,1,-1,1])