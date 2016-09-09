function [X,W]=Imp_quad_2D(y,lam,thres,xl,xu,method)

%% first is gaussian importance quadrature
% if strcmp(method,'gaussian1')==1
%     NN=15; % the number of points on x1 axis
%     ny=size(y,1);
%     x1=linspace(xl(1),xu(1),NN)';
%     X1=repmat(x1,1,ny);
%     Lam=repmat(lam',NN,1);
%     ex1=y(:,1)';
%     ex2=y(:,2)';
%     X2=Lam.*(X1.^repmat(ex1,NN,1));
%     uex2=unique(ex2);
%     uex2=uex2(end:-1:1);
%     Y2=zeros(NN,length(uex2));
%     for i=1:1:length(uex2)
%         ind=find(uex2(i)==ex2);
%         Y2(:,i)=sum(X2(:,ind),2);
%         if uex2(i)==0
%             Y2(:,i)=Y2(:,i)-log(thres);
%         end
%     end
%     P=[];
%     for i=1:1:NN
%     p=roots(Y2(i,:));
%     p(find(imag(p)~=0))=[];
%     if isempty(p)==0
%     P=vertcat(P,[repmat(x1(i),length(p),1),p]);
% %     P=vertcat(P,[x1(i),max(p)]);
%     end
%     end
%     
%% take gauss comp of max ent pdf
    if 1
        A=-2*real([lam(4),lam(5)/2;lam(5)/2,lam(6)]);
        b=[lam(2),lam(3)]';
        mu=inv(A)*b
        if isnan(sum(mu))
            mu=[0;0];
            P=eye(2);
            A=eye(2);
        else
        P=0.5*inv(A)
        end
        
        C=lam(1)+0.5*b'*inv(A)*b;
        NN=@(x)1/sqrt(det(2*pi*P))*exp(-1/2*(x-mu)'*inv(P)*(x-mu));
        [X,W]= GH_points(mu,P,20);
        for i=1:1:length(W)
        W(i)=W(i)/NN(X(i,:)');
        end
    end
    
end