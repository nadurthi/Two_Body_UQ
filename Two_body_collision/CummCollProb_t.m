function f=CummCollProb_t(R,mur,muv,Pr,Prv,Pv,YY,meth,X1vel,X2vel,w1,w2)
% YY is the max-ent pdf
mur=mur(:);
muv=muv(:);

[leb] = getLebedevSphere(590);
Xsph=[leb.x,leb.y,leb.z];
wsph=leb.w;

v0nc=@(nc)nc(:)'*(muv+Prv'*inv(Pr)*(R*nc(:)-mur));
sig2nc=@(nc)nc(:)'*(Pv-Prv'*inv(Pr)*Prv)*nc(:);

vnc=@(nc)sqrt(sig2nc(nc))/sqrt(2*pi)*exp(-v0nc(nc)^2/(2*sig2nc(nc)))-(v0nc(nc)/2)*(1-erf(v0nc(nc)/(sqrt(sig2nc(nc))*sqrt(2))));

% Now cal averaegt vel with quad pts without using the conditional pdf
n=size(X1vel,2);
NN=length(w1)*length(w2);
W=zeros(NN,1);
Xd=zeros(NN,n);
k=1;
for i=1:1:length(w1)
    for j=1:1:length(w2)
        Xd(k,:) =X1vel(i,:)-X2vel(j,:);
        W(k)=w1(i)*w2(j);
        k=k+1;
    end
end



%===============================================


Y=YY{1};
lam=YY{2};
muu=YY{3};
iP=YY{4};
% pr=YY{5};
muu=muu(:)';
pr=@(X)pdf_MaxEnt((iP*(X(:)'-muu)')',lam,Y)*det(iP);


f=0;
for i=1:1:length(wsph)
    nc=Xsph(i,:)';
    
    
    F=sum(Xd.*repmat(nc',NN,1),2);
    ind=find(F<=0);
    F=abs(sum(Xd(ind,:).*repmat(nc',length(ind),1),2));
    velavg=sum(F.*W(ind));
%     velavg=0;
%     for j=1:1:NN
%         vi=Xd(j,:)';
%         if dot(vi,nc)<=0
%         velavg=velavg+W(j)*abs(dot(vi,nc));
%         end
%     end
    
    
    if strcmp(meth,'gauss')
        pp=mvnpdf(R*nc,mur,Pr);
        f=f+wsph(i)*R^2*pp*(1*vnc(nc));
    elseif strcmp(meth,'nongauss')
%         f=f+wsph(i)*R^2*pr(R*nc)*(1*vnc(nc));
        f=f+wsph(i)*R^2*pr(R*nc)*velavg;
    end
    
end


end