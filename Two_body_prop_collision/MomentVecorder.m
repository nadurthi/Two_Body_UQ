function ind=MomentVecorder(ee)
ee=ee(:)';
Y1=[     1     0     0
        0     1     0
        0     0     1];
% 
Y2=  [2     0     0
     1     0     1
     1     1     0
     0     1     1
     0     0     2
     0     2     0];
%  
Y3=[ 3     0     0
     2     0     1
     2     1     0
     1     1     1
     1     0     2
     1     2     0
     0     2     1
     0     1     2
     0     0     3
     0     3     0];
% 
% 
Y4=[ 4     0     0
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
% 
% 
Y5=[ 5     0     0 
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
% 
Y6=[6     0     0
     5     0     1
     5     1     0
     4     1     1
     4     0     2
     4     2     0
     3     2     1
     3     1     2
     3     0     3
     3     3     0
     2     3     1
     2     2     2
     2     1     3
     2     0     4
     2     4     0
     1     4     1
     1     3     2
     1     2     3
     1     1     4
     1     0     5
     1     5     0
     0     5     1
     0     4     2
     0     3     3
     0     2     4
     0     1     5
     0     0     6
     0     6     0];

ind=find(sum(abs(Y1-repmat(ee,size(Y1,1),1)),2)==0);
if length(ind)>0
    ind=[1,ind(1)];
    return;
end

Y=Y2;
ind=find(sum(abs(Y-repmat(ee,size(Y,1),1)),2)==0);
if length(ind)>0
    ind=[2,ind(1)];
    return;
end

Y=Y3;
ind=find(sum(abs(Y-repmat(ee,size(Y,1),1)),2)==0);
if length(ind)>0
    ind=[3,ind(1)];
    return;
end

Y=Y4;
ind=find(sum(abs(Y-repmat(ee,size(Y,1),1)),2)==0);
if length(ind)>0
    ind=[4,ind(1)];
    return;
end

Y=Y5;
ind=find(sum(abs(Y-repmat(ee,size(Y,1),1)),2)==0);
if length(ind)>0
    ind=[5,ind(1)];
    return;
end

Y=Y6;
ind=find(sum(abs(Y-repmat(ee,size(Y,1),1)),2)==0);
if length(ind)>0
    ind=[6,ind(1)];
    return;
end