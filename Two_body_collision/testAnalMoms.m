X11cut8=zeros(size(X1cut8,3),3);
X22cut8=zeros(size(X1cut8,3),3); 

for j=1:1:size(X1cut8,3)
        X11cut8(j,:)=X1cut8(50,1:3,j);
        X22cut8(j,:)=X2cut8(50,1:3,j);
end
    
D=zeros(size(X11cut8,1)*size(X22cut8,1),3);
 W=zeros(size(D,1),1);
 k=1;
 for i=1:1:size(X11cut8,1)
     for j=1:1:size(X22cut8,1)
     D(k,:)=X11cut8(i,:)-X22cut8(j,:);
     W(k)=w1cut8(i)*w2cut8(j);
     k=k+1;
     end
 end

[y1,M1x]=Cal_moments_samples(D ,W,1,'raw');
[y2,M2x]=Cal_moments_samples(D ,W,2,'raw');
[y3,M3x]=Cal_moments_samples(D ,W,3,'raw');
[y4,M4x]=Cal_moments_samples(D ,W,4,'raw');
[y5,M5x]=Cal_moments_samples(D ,W,5,'raw');
[y6,M6x]=Cal_moments_samples(D ,W,6,'raw');
 MMx.M1=M1x;
 MMx.M2=M2x;
 MMx.M3=M3x;
 MMx.M4=M4x;
 MMx.M5=M5x;
 MMx.M6=M6x;

[muu,P]=ptswts2muP(D,W);

Dt=zeros(size(D));
for i=1:1:size(D,1)
   Dt(i,:)=sqrtm(inv(P))*(D(i,:)-muu)'; 
end
[y1,M1tDth]=Cal_moments_samples(Dt ,W,1,'raw');
 [y2,M2tDth]=Cal_moments_samples(Dt ,W,2,'raw');
 [y3,M3tDth]=Cal_moments_samples(Dt ,W,3,'raw');
 [y4,M4tDth]=Cal_moments_samples(Dt ,W,4,'raw');
  [y5,M5tDth]=Cal_moments_samples(Dt ,W,5,'raw');
[y6,M6tDth]=Cal_moments_samples(Dt ,W,6,'raw');
 MDt={M1tDth';M2tDth';M3tDth';M4tDth';M5tDth';M6tDth'};

 
 
muu=Dcut8{1}.M1;
P=[Dcut8{1}.M2(1),Dcut8{1}.M2(3),Dcut8{1}.M2(2);Dcut8{1}.M2(3),Dcut8{1}.M2(6),Dcut8{1}.M2(4);Dcut8{1}.M2(2),Dcut8{1}.M2(4),Dcut8{1}.M2(5)];
P=P-muu'*muu;
MDtccr=AffineTransMoms(Dcut8{1},6,-sqrtm(inv(P))*muu',sqrtm(inv(P)));

[MDtccr.M1', MDt{1}']
[MDtccr.M2'- MDt{2}']
[MDtccr.M3'- MDt{3}']
[MDtccr.M4'- MDt{4}']
[MDtccr.M5'- MDt{5}']
[MDtccr.M6'- MDt{6}']



%% 
Dcut8{1}.M1'-M1x
Dcut8{1}.M2'-M2x
Dcut8{1}.M3'-M3x
Dcut8{1}.M4'-M4x
%%

MMy=AffineTransMoms(MMx,6,-sqrtm(inv(P))*muu',sqrtm(inv(P)));

[MMy.M1'-MDt{1}']
[MMy.M2'-MDt{2}']
[MMy.M3'-MDt{3}']
[MMy.M4'-MDt{4}']
[MMy.M5'-MDt{5}']
[MMy.M6'-MDt{6}']



%%
A=randn(6,6);
P1=(A*A').^2*1;
mu1=[10,80,60,10,10,10]';
% [X11cut8,w1cut8]=conjugate_dir_gausspts_till_8moment([10,10,10]',P1);
[X11cut,w1cut]=GH_points(mu1,P1,3);
A=randn(6,6);
P2=(A*A').^2*1;
mu2=[5,83,91,10,-10,10]';
[X22cut,w2cut]=GH_points(mu2,P2,3);
% [X22cut8,w2cut8]=conjugate_dir_gausspts_till_8moment([10,-10,10]',P2);
X11cut(:,[4,5,6])=[];
X22cut(:,[4,5,6])=[];

SS=abs(mean([mean(X11cut8,1);mean(X22cut8,1)],1));
X11cut=X11cut8;%./repmat(SS,length(w1cut8),1);
X22cut=X22cut8;%./repmat(SS,length(w1cut8),1);
w1cut=w1cut8;
w2cut=w2cut8;
% X11cut=X11cut(:,[1,2,3]);
% X22cut=X22cut(:,[1,2,3]);


mudiff=mu1-mu2;
Pdiff=P1+P2;
[XD,wD]=GH_points(mudiff,Pdiff,5);
XD=XD(:,[1,2,3]);

D=zeros(size(X11cut,1)*size(X22cut,1),3);
 W=zeros(size(D,1),1);
 k=1;
 for i=1:1:size(X11cut,1)
     for j=1:1:size(X22cut,1)
     D(k,:)=X11cut(i,:)-X22cut(j,:);
     W(k)=w1cut(i)*w2cut(j);
     k=k+1;
     end
 end

[y1,M1x]=Cal_moments_samples(D ,W,1,'raw');
[y2,M2x]=Cal_moments_samples(D ,W,2,'raw');
[y3,M3x]=Cal_moments_samples(D ,W,3,'raw');
[y4,M4x]=Cal_moments_samples(D ,W,4,'raw');
[y5,M5x]=Cal_moments_samples(D ,W,5,'raw');
[y6,M6x]=Cal_moments_samples(D ,W,6,'raw');
 MMx.M1=M1x;
 MMx.M2=M2x;
 MMx.M3=M3x;
 MMx.M4=M4x;
 MMx.M5=M5x;
 MMx.M6=M6x;
 
[y1,M1x]=Cal_moments_samples(X11cut ,w1cut,1,'raw');
[y2,M2x]=Cal_moments_samples(X11cut ,w1cut,2,'raw');
[y3,M3x]=Cal_moments_samples(X11cut ,w1cut,3,'raw');
[y4,M4x]=Cal_moments_samples(X11cut ,w1cut,4,'raw');
[y5,M5x]=Cal_moments_samples(X11cut ,w1cut,5,'raw');
[y6,M6x]=Cal_moments_samples(X11cut,w1cut,6,'raw');
Mx.M1=M1x;
Mx.M2=M2x;
Mx.M3=M3x;
Mx.M4=M4x;
Mx.M5=M5x;
Mx.M6=M6x;


[y1,M1y]=Cal_moments_samples(X22cut ,w2cut,1,'raw');
[y2,M2y]=Cal_moments_samples(X22cut ,w2cut,2,'raw');
[y3,M3y]=Cal_moments_samples(X22cut ,w2cut,3,'raw');
[y4,M4y]=Cal_moments_samples(X22cut ,w2cut,4,'raw');
[y5,M5y]=Cal_moments_samples(X22cut ,w2cut,5,'raw');
[y6,M6y]=Cal_moments_samples(X22cut ,w2cut,6,'raw');
My.M1=M1y;
My.M2=M2y;
My.M3=M3y;
My.M4=M4y;
My.M5=M5y;
My.M6=M6y;


MD=MissMomsFromIndividMoms(Mx,My,6);

[MMx.M1,MD.M1',MDg3.M1]
[MMx.M2,MD.M2',MDg3.M2]
[MMx.M3,MD.M3',MDg3.M3]
[MMx.M4,MD.M4',MDg3.M4]
[MMx.M5,MD.M5',MDg3.M5]
[MMx.M6,MD.M6',MDg3.M6]

[MMx.M4-MDg3.M4,MD.M4'-MDg3.M4]
[MMx.M5-MDg3.M5,MD.M5'-MDg3.M5]
[MMx.M6-MDg3.M6,MD.M6'-MDg3.M6]

[y1,M1Dg]=Cal_moments_samples(XD,wD,1,'raw');
[y2,M2Dg]=Cal_moments_samples(XD,wD,2,'raw');
[y3,M3Dg]=Cal_moments_samples(XD,wD,3,'raw');
[y4,M4Dg]=Cal_moments_samples(XD,wD,4,'raw');
[y5,M5Dg]=Cal_moments_samples(XD,wD,5,'raw');
[y6,M6Dg]=Cal_moments_samples(XD,wD,6,'raw');
MDg.M1=M1Dg;
MDg.M2=M2Dg;
MDg.M3=M3Dg;
MDg.M4=M4Dg;
MDg.M5=M5Dg;
MDg.M6=M6Dg;

[y1,M1Dg]=Cal_moments_samples(XD,wD,1,'raw');
[y2,M2Dg]=Cal_moments_samples(XD,wD,2,'raw');
[y3,M3Dg]=Cal_moments_samples(XD,wD,3,'raw');
[y4,M4Dg]=Cal_moments_samples(XD,wD,4,'raw');
[y5,M5Dg]=Cal_moments_samples(XD,wD,5,'raw');
[y6,M6Dg]=Cal_moments_samples(XD,wD,6,'raw');
MDg3.M1=M1Dg;
MDg3.M2=M2Dg;
MDg3.M3=M3Dg;
MDg3.M4=M4Dg;
MDg3.M5=M5Dg;
MDg3.M6=M6Dg;

[MDg.M1(1:3),MDg3.M1]
[MDg.M2(1:3),MDg3.M2]
[MDg.M3(1:3),MDg3.M3]
[MDg.M4(1:3),MDg3.M4]
[MDg.M5(1:3),MDg3.M5]
[MDg.M6(1:3),MDg3.M6]



MMy=AffineTransMoms(MMx,6,-sqrtm(inv(P))*muu',sqrtm(inv(P)));



[YJ,MJ]=AnalIndependentJointMoms6D(Mx,My,6);

[MD,MD6] =missdistmomsMODIFIED(Mx,My,6);

[MD.y1,MMx.M1]
[MD.y2,MMx.M2]
[MD.y3,MMx.M3]
[MD.y4-MMx.M4]
[MD.y5,MMx.M5]
[MD.y6,MMx.M6]

%%
syms x1 x2 x3 y1 y2 y3
d1=x1-y1;
d2=x2-y2;
d3=x3-y3;
d4=x1;
d5=x2;
d6=x3;

fid = fopen('AnalMissMoms6D.m', 'w');
Y1=YJ.y1;
for i=1:1:size(Y1,1)
fprintf(fid,'%s ;\n',strcat('M1(',num2str(i),')=',char(expand(d1^Y1(i,1)*d2^Y1(i,2)*d3^Y1(i,3)*d4^Y1(i,4)*d5^Y1(i,5)*d6^Y1(i,6)))))
end
Y2=YJ.y2;
for i=1:1:size(Y2,1)
fprintf(fid,'%s ;\n',strcat('M2(',num2str(i),')=',char(expand(d1^Y2(i,1)*d2^Y2(i,2)*d3^Y2(i,3)*d4^Y2(i,4)*d5^Y2(i,5)*d6^Y2(i,6)))))
end
Y3=YJ.y3;
for i=1:1:size(Y3,1)
fprintf(fid,'%s ;\n',strcat('M3(',num2str(i),')=',char(expand(d1^Y3(i,1)*d2^Y3(i,2)*d3^Y3(i,3)*d4^Y3(i,4)*d5^Y3(i,5)*d6^Y3(i,6)))))
end
Y4=YJ.y4;
for i=1:1:size(Y4,1)
fprintf(fid,'%s ;\n',strcat('M4(',num2str(i),')=',char(expand(d1^Y4(i,1)*d2^Y4(i,2)*d3^Y4(i,3)*d4^Y4(i,4)*d5^Y4(i,5)*d6^Y4(i,6)))))
end
Y5=YJ.y5;
for i=1:1:size(Y5,1)
fprintf(fid,'%s ;\n',strcat('M5(',num2str(i),')=',char(expand(d1^Y5(i,1)*d2^Y5(i,2)*d3^Y5(i,3)*d4^Y5(i,4)*d5^Y5(i,5)*d6^Y5(i,6)))))
end
Y6=YJ.y6;
for i=1:1:size(Y6,1)
fprintf(fid,'%s ;\n',strcat('M6(',num2str(i),')=',char(expand((d1^Y6(i,1))*(d2^Y6(i,2))*d3^Y6(i,3)*d4^Y6(i,4)*d5^Y6(i,5)*d6^Y6(i,6)))))
end
fclose(fid)

sum(w1cut.*prod(X11cut.^repmat([0,0,4],745,1),2))
Mx.M4(14)
p1=1.382578600940526e+015