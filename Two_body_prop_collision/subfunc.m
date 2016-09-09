
function Moms=subfunc(D,W,n)

y1=[ 1     0     0
     0     1     0
     0     0     1];

 y2=[2     0     0
     1     0     1
     1     1     0
     0     1     1
     0     0     2
     0     2     0];
y3=[3     0     0
     2     0     1
     2     1     0
     1     1     1
     1     0     2
     1     2     0
     0     2     1
     0     1     2
     0     0     3
     0     3     0];
 
 y4=[4     0     0
     3     0     1
     3     1     0
     2     1     1
     2     0     2
     2     2     0
     1     2     1
     1     1     2
     1     0     3
     1     3     0
     0     3     1
     0     2     2
     0     1     3
     0     0     4
     0     4     0];
 y5=[ 5     0     0 
     4     0     1
     4     1     0
     3     1     1
     3     0     2
     3     2     0
     2     2     1
     2     1     2
     2     0     3
     2     3     0
     1     3     1
     1     2     2
     1     1     3
     1     0     4
     1     4     0
     0     4     1
     0     3     2
     0     2     3
     0     1     4
     0     0     5
     0     5     0];
 
 Moms.M1=zeros(1,size(y1,1));
 Moms.M2=zeros(1,size(y2,1));
 Moms.M3=zeros(1,size(y3,1));
 Moms.M4=zeros(1,size(y4,1));
 Moms.M5=zeros(1,size(y5,1));
 
     for k=1:1:size(y1,1)
         Y=repmat(y1(k,:),n,1);
         Moms.M1(k)=Moms.M1(k)+sum(W.*prod(D.^Y,2));
     end
     for k=1:1:size(y2,1)
         Y=repmat(y2(k,:),n,1);
    Moms.M2(k)=Moms.M2(k)+sum(W.*prod(D.^Y,2));
     end
     for k=1:1:size(y3,1)
         Y=repmat(y3(k,:),n,1);
    Moms.M3(k)=Moms.M3(k)+sum(W.*prod(D.^Y,2));
     end
     for k=1:1:size(y4,1)
         Y=repmat(y4(k,:),n,1);
    Moms.M4(k)=Moms.M4(k)+sum(W.*prod(D.^Y,2));
     end
     for k=1:1:size(y5,1)
         Y=repmat(y5(k,:),n,1);
    Moms.M5(k)=Moms.M5(k)+sum(W.*prod(D.^Y,2));
     end
     
end