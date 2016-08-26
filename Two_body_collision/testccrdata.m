N=100000;
X1_1=zeros(N,3);
X1_2=zeros(N,3);
load('CollpaperSimsMC_2')
for j=1:1:N
X1_1(j,:)=X1mc(1,1:3,j);
X1_2(j,:)=X2mc(1,1:3,j);
end
X2_1=zeros(N,3);
X2_2=zeros(N,3);
load('CollpaperSimsMC_3')
for j=1:1:N
X2_1(j,:)=X1mc(1,1:3,j);
X2_2(j,:)=X2mc(1,1:3,j);
end
max(max(abs(X1_1-X2_1)))