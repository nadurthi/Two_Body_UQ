function H=plot_enttropy_sats(XsigSat,k,coll,Tvec)
nsat=size(XsigSat,1);
H=zeros(k,nsat);
s={};
for i=1:1:nsat
    for t=1:1:k
    P=XsigSat{i,2}(t,:);
    n=sqrt(length(P));
%     H(ks+1,i)=n/2*(1+log(2*pi))+0.5*log(det(reshape(P,n,n)));
    H(t,i)=norm(reshape(P,n,n),'fro');
    end
    s{i}=strcat('Sat ',num2str(i));
end
t=Tvec(1:k)/3600;
plot(repmat(t',1,nsat),H,'linewidth',2)
% plot(t',max(H,[],2),coll,'linewidth',2)
% legend('Sum of Frob norm')
legend(s,'Location','NorthWest')
end