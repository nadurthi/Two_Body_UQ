function ERR=plot_est_err_sats_CUMULATE_orb(XsigSat,k,coll,Tvec,ytruth_orb)
Re=6378.1;
nnn=length(XsigSat{1,1}(1,:));
nsat=size(XsigSat,1);
H=zeros(k,nsat);
s={};
err=zeros(k,nsat);
for i=1:1:nsat
   for t=1:1:k
%     err(t,i)=sqrt(1/t*sum(sqrt(1/nnn*sum((XsigSat{i,1}(1:t,:)-ytruth{i}(1:t,:)).^2)).^2));
   OE = XYZ2OE(XsigSat{i,1}(t,:));
   OE=OE(:)';
    err(t,i)=sqrt(sum((OE./[Re,1,1,1,1,1]-ytruth_orb{i}(t,:)./[Re,1,1,1,1,1]).^2)/nnn);
   end
end
t=Tvec(1:k)/3600;
errmax=max(err,[],2);
errmin=min(err,[],2);
errmean=mean(err,2);
ERR=[errmax,errmin,errmean];
% figure
% plot(t,errmin,'r--',t,errmax,'b--',t,errmean,'k','linewidth',2)
% figure
% plot(t,err)

