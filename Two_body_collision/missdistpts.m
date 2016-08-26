function MoMs=missdistpts(Xx,wx,Xy,wy)


 
MOMS=cell(size(Xx,1),1);

 for i=1:size(Xx,1)
     %for j=1:1:size(Xy,1)
     n=size(Xy,1);
     D=repmat(Xx(i,:),n,1)-Xy;
     W=repmat(wx(i),n,1).*wy;
         
      MOMS{i}=subfunc(D,W,n);
     %end
 end
 
 MoMs.M1=zeros(1,3);
 MoMs.M2=zeros(1,6);
 MoMs.M3=zeros(1,10);
 MoMs.M4=zeros(1,15);
 MoMs.M5=zeros(1,21);
 
 for i=1:1:size(Xx,1)
    MoMs.M1=MoMs.M1+MOMS{i}.M1;
    MoMs.M2=MoMs.M2+MOMS{i}.M2;
    MoMs.M3=MoMs.M3+MOMS{i}.M3;
    MoMs.M4=MoMs.M4+MOMS{i}.M4;
    MoMs.M5=MoMs.M5+MOMS{i}.M5;
end
