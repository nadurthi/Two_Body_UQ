function H=plot_enttropy_sats_CUMULATE(XsigSat,k,coll,Tvec,further)
nsat=size(XsigSat,1);
H=zeros(k,nsat);
s={};
for i=1:1:nsat
    for t=1:1:k
    P=XsigSat{i,2}(t,:);
    n=sqrt(length(P));
%     H(ks+1,i)=n/2*(1+log(2*pi))+0.5*log(det(reshape(P,n,n)));
    H(t,i)=norm(reshape(P,n,n),'fro');
    H(t,i)=sqrt(max(eig(reshape(P,n,n))));
    end
    s{i}=strcat('Sat ',num2str(i));
end
t=Tvec(1:k)/3600;
Hmax=max(H,[],2);
Hmin=min(H,[],2);
Hmean=mean(H,2);
HH=[Hmax,Hmin,Hmean];
if max(strcmpi(further,'plotit'))==1
plot(t,Hmin,'r--',t,Hmax,'b--',t,Hmean,'k','linewidth',2)
legend('min','max','mean')
end

end