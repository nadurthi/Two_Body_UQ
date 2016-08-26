function [M1,M2,M3,M4]=ORBmomsMC(ss)
mu=398601.2;
load(strcat('eg1CollisionBodies',num2str(ss)))
load('yMOMS')

M1=zeros(110,size(yMOMS{1},1));
M2=zeros(110,size(yMOMS{2},1));
M3=zeros(110,size(yMOMS{3},1));
M4=zeros(110,size(yMOMS{4},1));
X=zeros(3,1);
V=zeros(3,1);

for t=1:1:110
    for i=1:1:size(X1mc,3)
        X(1)=X1mc(t,1,i)';
        X(2)=X1mc(t,2,i)';
        X(3)=X1mc(t,3,i)';
        
        V(1)=X1mc(t,4,i)';
        V(2)=X1mc(t,5,i)';
        V(3)=X1mc(t,6,i)';
        [a,e,E,w,inc,Om] = XYZ2OE(X,V,mu);
        
        for j=1:1:size(yMOMS{1},1)
            M1(t,j)=M1(t,j)+1/(50*100000)*prod([a,e,E,w,inc,Om].^yMOMS{1}(j,:));
        end
        for j=1:1:size(yMOMS{2},1)
            M2(t,j)=M2(t,j)+1/(50*100000)*prod([a,e,E,w,inc,Om].^yMOMS{2}(j,:));
        end
        for j=1:1:size(yMOMS{3},1)
            M3(t,j)=M3(t,j)+1/(50*100000)*prod([a,e,E,w,inc,Om].^yMOMS{3}(j,:));
        end
        for j=1:1:size(yMOMS{4},1)
            M4(t,j)=M4(t,j)+1/(50*100000)*prod([a,e,E,w,inc,Om].^yMOMS{4}(j,:));
        end
        
        
    end
end