DVcut8=cell(31,1);
t=40:70;
for i=t
    X11cut8=zeros(size(X1cut8,3),6);
    X22cut8=zeros(size(X1cut8,3),6);

    for j=1:1:size(X1cut8,3)
        X11cut8(j,:)=X1cut8(i,1:6,j);
        X22cut8(j,:)=X2cut8(i,1:6,j);
    end
    
    
    DVcut8{i-39}=missdistpts(X11cut8,w1cut8,X22cut8,w2cut8);
    
end
MV1=zeros(length(t),length(DVcut8{1}.M1));
MV2=zeros(length(t),length(DVcut8{1}.M2));
MV3=zeros(length(t),length(DVcut8{1}.M3));
MV4=zeros(length(t),length(DVcut8{1}.M4));

for k=1:1:length(t)
MV1(k,:)=DVcut8{k}.M1;
MV2(k,:)=DVcut8{k}.M2;
MV3(k,:)=DVcut8{k}.M3;
MV4(k,:)=DVcut8{k}.M4;
end



i=1
load('CollpaperSimsMC_20')
X11mc=zeros(size(X1mc,3),6);
for j=1:1:size(X1mc,3)
X11mc(j,:)=X1mc(i,1:6,j);
end


load('CollpaperSimsMC_12')
X22mc=zeros(size(X1mc,3),6);
for j=1:1:size(X1mc,3)
X22mc(j,:)=X1mc(i,1:6,j);
end

max(max(abs(X11mc-X22mc)))