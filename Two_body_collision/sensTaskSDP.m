% clc
clear
close all
% P=magic(5)*magic(5)'/100;
Px=[2,-1;-1,4.5];
mx=[3,5];
Z1=@(x)[norm([x(1)-2,x(2)-2]);atan2(x(2)-2,x(1)-2)];
Z2=@(x)[norm([x(1)+2,x(2)+2]);atan2(x(2)+2,x(1)+2)];
Z3=@(x)[norm([x(1)-10,x(2)]);atan2(x(2),x(1)-10)];
hjack=@(x,a,b)[(x(1)-a)/norm([x(1)-a,x(2)-b]),(x(2)-b)/norm([x(1)-a,x(2)-b]);(b-x(2))/norm([x(1)-a,x(2)-b])^2,(-a+x(1))/norm([x(1)-a,x(2)-b])^2];
R1=diag([1^2,(1*pi/180)^2]);
R2=diag([5^2,(0.1*pi/180)^2]);
R3=diag([0.5^2,(10*pi/180)^2]);
hn=2;
nz=3;
H=[hjack(mx,2,2);hjack(mx,-2,-2);hjack(mx,10,0)];
H1=hjack(mx,2,2);
H2=hjack(mx,-2,-2);
H3=hjack(mx,10,0);



% [xut,wut]=UT_sigmapoints(mx',Px,2);
% z1ut=zeros(length(wut),2);
% for i=1:1:length(wut)
%     z1ut(i,:)=Z1(xut(i,:));
% end

[x,w]=GH_points(mx',Px,5);
z1=zeros(length(w),2);
z2=zeros(length(w),2);
z3=zeros(length(w),2);
for i=1:1:length(w)
    z1(i,:)=Z1(x(i,:));
    z2(i,:)=Z2(x(i,:));
    z3(i,:)=Z3(x(i,:));
end

Z=[z1,z2,z3];
Hlin=StatLinearize_sigmapts(Z,x,w);
H1lin=Hlin(1:2,:);
H2lin=Hlin(3:4,:);
H3lin=Hlin(5:6,:);

[mz,Pz]=MeanCov(Z,w);
Pz=Pz+blkdiag(R1,R2,R3);
Pxz=CrossCov(x,mx,Z,mz,w);


v=zeros(3,1);
v(1)=1;
v(2)=1;
v(3)=0;


% True infor
if v(1)==1 && v(2)==1 && v(3)==1
    InfoT=Px-Pxz*inv(Pz)*Pxz';
elseif v(1)==0 && v(2)==0 && v(3)==0
    InfoT=Px;
else
    
    R=cell(3,1);
    R{1}=R1;R{2}=R2;R{3}=R3;
    RR=[];
    Zind=[];
    for i=1:1:3
        if v(i)==1
            Zind=horzcat(Zind,Z(:,(i-1)*hn+1:i*hn));
            RR=blkdiag(RR,R{i});
            
        end
    end
    %     keyboard
    [mmz,PPz]=MeanCov(Zind,w);
    PPz=PPz+RR;
    PPxz=CrossCov(x,mx,Zind,mmz,w);
    
    InfoT=0.5*log(det(Px)/det(Px-PPxz*inv(PPz)*PPxz'));
end
%FIM optimization
F(1)=trace(H1lin'*inv(R1)*H1lin);
F(2)=trace(H2lin'*inv(R2)*H2lin);
F(3)=trace(H3lin'*inv(R3)*H3lin);
vfim=bintprog(-F',[],[],[1 1 1],2)

F(1)=trace(H1'*inv(R1)*H1);
F(2)=trace(H2'*inv(R2)*H2);
F(3)=trace(H3'*inv(R3)*H3);
vfim=bintprog(-F',[],[],[1 1 1],2)


[log(det(Px)/det(InfoT));trace(InfoT)]
% [det(Info),det(InfoT);trace(Info),trace(InfoT);]
%%
[110]
1.8968
22.388

[101]
3.3299
13.143

[011]
2.8906
16.546

[111]
3.7373
11.536

[100]
1.2834
29.947


[010]
1.6475
25.026
[001]
2.2696
23.516
%%

% Then we define the optimization problem
mpol('v',1,3)
[p,q]=MI_num_dem_poly_zver(v,Px,Pxz,Pz,hn,nz,0);
K=(eye(3)-diag(v))*diag(v);
K=diag(K);
K=[K==0;sum(v)==1];
Prob = msdp(min(q),K);
[status,obj] = msol(Prob) ;
vsolz = round(double(v))

% pause % Strike any key to continue.
[p,q]=MI_num_dem_value([1 0 1],Px,Pxz,Pz,hn,nz,0)
mpol('v',1,3)
[p,q]=MI_num_dem_poly(v,Px,Pxz,Pz,hn,nz,0);
K=(eye(3)-diag(v))*diag(v);
K=diag(K);
K=[K==0;sum(v)==1];
Prob = msdp(min(p),mom(q)==1,K);
[status,obj] = msol(Prob) ;
vsolx = double(v)


cvx_begin sdp
variable x(3) nonnegative
maximize(log_det(inv(Px)+x(1)*H1lin'*inv(R1)*H1lin+x(2)*H2lin'*inv(R2)*H2lin+x(3)*H3lin'*inv(R3)*H3lin))
subject to
sum(x)==2
x(1)<=1
x(2)<=1
x(3)<=1
cvx_end
x

cvx_begin sdp
variable y(3)
variable x(2)
maximize(log_det(inv(Px)+x(1)*H1'*inv(R1)*H1+x(2)*H2'*inv(R2)*H2+x(3)*H3'*inv(R3)*H3))
subject to
[x(1),x(2),1;1,2,1;-3,5,1]*y==[]
cvx_end
x


y=coef(qd);
[ys,ind]=sort(y);
cvx_begin
variable x(n) nonnegative
minimize(norm(ys-x,1)+norm(A*x,4))
subject to
ys>=x
cvx_end
for i=1:1:7
   yr(i)=x(ind==i) ;
end
[y,yr']
%%
% plotting
[xx,yy]=meshgrid(-30:0.1:30);
pp=zeros(size(xx));

for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx,Px);
    end
end
figure
contour(xx,yy,pp,30);
hold on
plot(mx(1),mx(2),'ro',2,2,'ks',-2,-2,'ks',10,0,'ks')
plot(x(:,1),x(:,2),'b+')


%%
MCx=mvnrnd(mx,Px,10000);
MCz=zeros(size(MCx,1),2);
for i=1:1:size(MCx,1)
    MCz(i,:)=Z1(MCx(i,:));
end
mean(MCz,1)


%%
%    Multiple targets
Px1=[20,-10;-10,45];
mx1=[3,5];
Px2=[30,15;15,10];
mx2=[-3,5];
a1=2;
b1=2;
a2=-2;
b2=-2;
a3=10;
b3=0;
Z1=@(x)[norm([x(1)-a1,x(2)-b1]);atan2(x(2)-b1,x(1)-a1)];
Z2=@(x)[norm([x(1)-a2,x(2)-b2]);atan2(x(2)-b2,x(1)-a2)];
Z3=@(x)[norm([x(1)-a3,x(2)-b3]);atan2(x(2)-b3,x(1)-a3)];
hjack=@(x,a,b)[(x(1)-a)/norm([x(1)-a,x(2)-b]),(x(2)-b)/norm([x(1)-a,x(2)-b]);(b-x(2))/norm([x(1)-a,x(2)-b])^2,(-a+x(1))/norm([x(1)-a,x(2)-b])^2];
R1=diag([1^2,(1*pi/180)^2]);
R2=diag([5^2,(0.1*pi/180)^2]);
R3=diag([0.5^2,(10*pi/180)^2]);
hn=2;
fn=2;
nz=3;
nx=2;

[x1,w1]=GH_points(mx1',Px1,5);
[x2,w2]=GH_points(mx2',Px2,5);
z1_x1=zeros(length(w1),2);
z1_x2=zeros(length(w1),2);
z2_x1=zeros(length(w1),2);
z2_x2=zeros(length(w1),2);
z3_x1=zeros(length(w1),2);
z3_x2=zeros(length(w1),2);

for i=1:1:length(w1)
    z1_x1(i,:)=Z1(x1(i,:));
    z1_x2(i,:)=Z1(x2(i,:));
    
    z2_x1(i,:)=Z2(x1(i,:));
    z2_x2(i,:)=Z2(x2(i,:));
    
    z3_x1(i,:)=Z3(x1(i,:));
    z3_x2(i,:)=Z3(x2(i,:));
end

Z1=[z1_x1,z2_x1,z3_x1];
Z2=[z1_x2,z2_x2,z3_x2];

Hlin_x1=StatLinearize_sigmapts(Z1,x1,w1);
Hlin_x2=StatLinearize_sigmapts(Z2,x2,w2);

H1lin_x1=Hlin_x1(1:2,:);
H2lin_x1=Hlin_x1(3:4,:);
H3lin_x1=Hlin_x1(5:6,:);

H1lin_x2=Hlin_x2(1:2,:);
H2lin_x2=Hlin_x2(3:4,:);
H3lin_x2=Hlin_x2(5:6,:);

[mz1,Pz1]=MeanCov(Z1,w1);
Pz1=Pz1+blkdiag(R1,R2,R3);
Pxz1=CrossCov(x1,mx1,Z1,mz1,w1);

[mz2,Pz2]=MeanCov(Z2,w2);
Pz2=Pz2+blkdiag(R1,R2,R3);
Pxz2=CrossCov(x2,mx2,Z2,mz2,w2);

%% FIM optimization
H1_x1=hjack(mx1,a1,b1);
H2_x1=hjack(mx1,a2,b2);
H3_x1=hjack(mx1,a3,b3);

H1_x2=hjack(mx2,a1,b1);
H2_x2=hjack(mx2,a2,b2);
H3_x2=hjack(mx2,a3,b3);

% F(1,1)=trace(H1lin_x1'*inv(R1)*H1lin_x1);
% F(1,2)=trace(H2lin_x1'*inv(R2)*H2lin_x1);
% F(1,3)=trace(H3lin_x1'*inv(R3)*H3lin_x1);
%
% F(2,1)=trace(H1lin_x2'*inv(R1)*H1lin_x2);
% F(2,2)=trace(H2lin_x2'*inv(R2)*H2lin_x2);
% F(2,3)=trace(H3lin_x2'*inv(R3)*H3lin_x2);

F(1,1)=trace(H1_x1'*inv(R1)*H1_x1);
F(1,2)=trace(H2_x1'*inv(R2)*H2_x1);
F(1,3)=trace(H3_x1'*inv(R3)*H3_x1);

F(2,1)=trace(H1_x2'*inv(R1)*H1_x2);
F(2,2)=trace(H2_x2'*inv(R2)*H2_x2);
F(2,3)=trace(H3_x2'*inv(R3)*H3_x2);

vv=bintprog(-[F(1,:)';F(2,:)'],[],[],[eye(3),eye(3)],ones(3,1));
vfim=[vv(1:3)';vv(4:6)']

cvx_begin
variable mmu(2) binary
variable u(2,3) binary
maximize(mmu(1)*trace(inv(Px1))+u(1,1)*F(1,1)+u(1,2)*F(1,2)+u(1,3)*F(1,3)...
        +mmu(2)*trace(inv(Px2))+u(2,1)*F(2,1)+u(2,2)*F(2,2)+u(2,3)*F(2,3) )
subject to
 u(1,1)+u(2,1)==1
 u(1,2)+u(2,2)==1
 u(1,3)+u(2,3)==1
 1/3*(u(1,1)+u(1,2)+u(1,3))<=mmu(1)<=(u(1,1)+u(1,2)+u(1,3))
 1/3*(u(2,1)+u(2,2)+u(2,3))<=mmu(2)<=(u(2,1)+u(2,2)+u(2,3))
cvx_end

%%  Assume MI is sum of individual i.e. maximixe the bound
% eg I(x1:z1,z2,z3)=I(x1:z1)+I(x1:z2)+I(x1:z3)
tic
MI=zeros(2,3);

ZZ1=[z1_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT11=Px1-PPxz1*inv(PPz1)*PPxz1';
MI(1,1)=log(det(Px1)/det(InfoT11));

ZZ1=[z2_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R2);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT12=Px1-PPxz1*inv(PPz1)*PPxz1';
MI(1,2)=log(det(Px1)/det(InfoT12));

ZZ1=[z3_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT13=Px1-PPxz1*inv(PPz1)*PPxz1';
MI(1,3)=log(det(Px1)/det(InfoT13));

ZZ2=[z1_x2]
[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R1);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);
InfoT21=Px2-PPxz2*inv(PPz2)*PPxz2';
MI(2,1)=log(det(Px2)/det(InfoT21));

ZZ2=[z2_x2]
[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R2);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);
InfoT22=Px2-PPxz2*inv(PPz2)*PPxz2';
MI(2,2)=log(det(Px2)/det(InfoT22));

ZZ2=[z3_x2]
[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R3);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);
InfoT23=Px2-PPxz2*inv(PPz2)*PPxz2';
MI(2,3)=log(det(Px2)/det(InfoT23));

vv=bintprog(-[MI(1,:)';MI(2,:)'],[],[],[eye(3),eye(3)],ones(3,1));
vMI=[vv(1:3)';vv(4:6)']
toc
%% x1 z1 z2 z3
ZZ1=[z1_x1,z2_x1,z3_x1];

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1,R2,R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2;
clc
log(det(Px1)/det(InfoT1))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)
  3.7373
   39.0318
   [trace(InfoT1),trace(Px2)]
      [det(InfoT1),det(Px2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       19.055       11.536          800           65
%           800           65          800           65]

%% x2 z1 z2 z3
ZZ2=[z1_x2,z2_x2,z3_x2];

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R1,R2,R3);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1;
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px2)/det(InfoT2))
trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
3.3676
21.4029
   [trace(Px1),trace(InfoT2)]
      [det(Px1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [800           65       19.975       13.592
%  800           65          800           65]

%% x1 z1 z2 , x2 z3
ZZ1=[z1_x1,z2_x1];
ZZ2=[z3_x2];

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1,R2);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);


[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+blkdiag(R3);
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
4.2316
31.5397
   [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [120.04       22.388       45.496       25.119
%  800           65          800           65]
%% x1 z1 z3 , x2 z2
ZZ1=[z1_x1,z3_x1];
ZZ2=[z2_x2];
RR1=blkdiag(R1,R3);
RR2=blkdiag(R2);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
 3.7728
 28.0290
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [28.637       13.143       261.65       33.369
%  800           65          800           65]
%% x1 z2 z3 , x2 z1
ZZ1=[z2_x1,z3_x1];
ZZ2=[z1_x2];
RR1=blkdiag(R2,R3);
RR2=blkdiag(R1);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
3.6991
39.8977
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       44.436       16.546       132.02       30.473
%           800           65          800           65]

%% x1 z1 , x2 z2 Z3
ZZ1=[z1_x1];
ZZ2=[z2_x2,z3_x2];
RR1=blkdiag(R1);
RR2=blkdiag(R2,R3);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
4.887
20.5370
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       221.67       29.947       33.223       19.328
%           800           65          800           65]

%% x1 z2 , x2 z1 Z3
ZZ1=[z2_x1];
ZZ2=[z1_x2,z3_x2];
RR1=blkdiag(R2);
RR2=blkdiag(R1,R3);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';
clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
4.8962
32.4056
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       154.03       25.026       21.825       14.565
%           800           65          800           65]

%% x1 z3 , x2 z1 Z2

ZZ1=[z3_x1];
ZZ2=[z1_x2,z2_x2];
RR1=blkdiag(R3);
RR2=blkdiag(R1,R2);

[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+RR1;
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);

[mmz2,PPz2]=MeanCov(ZZ2,w2);
PPz2=PPz2+RR2;
PPxz2=CrossCov(x2,mx2,ZZ2,mmz2,w2);

InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';
InfoT2=Px2-PPxz2*inv(PPz2)*PPxz2';

clc
log(det(Px1)/det(InfoT1))+log(det(Px2)/det(InfoT2))
trace(PPz1)-trace(PPz1-PPxz1'*inv(Px1)*PPxz1)+trace(PPz2)-trace(PPz2-PPxz2'*inv(Px2)*PPxz2)
3.5786
28.8949
    [trace(InfoT1),trace(InfoT2)]
      [det(InfoT1),det(InfoT2)]
% [det(InfoT1),trace(InfoT1),det(InfoT2),trace(InfoT2);det(Px1),trace(Px1),det(Px2),trace(Px2)]
% [       82.686       23.516        101.5       25.149
%           800           65          800           65]


%%

Px=cell(2,1);
Pxz=cell(2,1);
Px{1}=Px1;  Px{2}=Px2;
Pxz{1}=Pxz1;  Pxz{2}=Pxz2;
Pz=cell(2,1);
Pz{1}=Pz1;  Pz{2}=Pz2;
mpol('v',nx,nz)
[p,q]=MI_num_dem_poly_zver_multi_x(v,Px,Pxz,Pz,hn,nz,nx,0)

K=(ones(nx,nz)-v).*v;
K=reshape(K,nx*nz,1);
mpol('S',nz,1);
for i=1:1:nz
    S(i)=sum(v(:,i));
end
K=[K==0;S==1];

Prob = msdp(max(q),K);
tic

[status,obj] = msol(Prob)
toc
vsolz = round(double(v))

%% is MI a sum of individual MI

% ZZ1=[z1_x1,z2_x1,z3_x1];3.7373
ZZ1=[z1_x1,z2_x1,z3_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1,R2,R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT1=Px1-PPxz1*inv(PPz1)*PPxz1';

ZZ1=[z1_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R1);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT11=Px1-PPxz1*inv(PPz1)*PPxz1';

ZZ1=[z2_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);

PPz1=PPz1+blkdiag(R2);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT12=Px1-PPxz1*inv(PPz1)*PPxz1';

ZZ1=[z3_x1]
[mmz1,PPz1]=MeanCov(ZZ1,w1);
PPz1=PPz1+blkdiag(R3);
PPxz1=CrossCov(x1,mx1,ZZ1,mmz1,w1);
InfoT13=Px1-PPxz1*inv(PPz1)*PPxz1';


log(det(Px1)/det(InfoT1))
log(det(Px1)/det(InfoT11))+log(det(Px1)/det(InfoT12))+log(det(Px1)/det(InfoT13))
%% matching the z-ver trace and det
% v=sym('v',[nx,nz]);
PPxz=Pxz2;
PPz=Pz2;
Px=Px2;
v=[1 0 0];

for i=1:1:nz
    PPxz(:,(i-1)*hn+1:i*hn)=v(i)*PPxz(:,(i-1)*hn+1:i*hn);
    for j=1:1:nz
%         if j~=i
            PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)=v(i)*PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn);
            PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn)=v(i)*PPz((j-1)*hn+1:j*hn,(i-1)*hn+1:i*hn);
%         end
    end
end

q=trace(PPz-PPxz'*inv(Px)*PPxz);
p=trace(PPz);

PPxz=Pxz2;
PPz=Pz2;
Px=Px2;
PPxz_red=[];
PPz_red=[];
for i=1:1:nz
    if v(i)==1
        PPxz_red=[PPxz_red,PPxz(:,(i-1)*hn+1:i*hn)];
    end
    PP=[];
    for j=1:1:nz
            if v(i)==1 && v(j)==1
                PP=[PP,PPz((i-1)*hn+1:i*hn,(j-1)*hn+1:j*hn)];
            end
    end
    PPz_red=[PPz_red;PP];
end
q_red=trace(PPz_red-PPxz_red'*inv(Px)*PPxz_red);
p_red=trace(PPz_red);

[q_red,q;
    p_red,p]

%% Plotting the targets
[xx,yy]=meshgrid(-20:0.1:20);
pp1=zeros(size(xx));
pp2=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp1(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx1,Px1);
        pp2(i,j)=mvnpdf([xx(i,j),yy(i,j)],mx2,Px2);
     end
end
contour(xx,yy,pp1)
hold on
contour(xx,yy,pp2)
plot(2,2,'ks',-2,-2,'ks',10,0,'ks','linewidth',4,'MarkerSize',8)

%% Plotting the targets
Px11_MI4= [5.7102    4.1758
    4.1758   10.8358];
Px22_MI4= [ 21.1347   11.0866
   11.0866    7.3967];

Px11_MUB=[11.4930    2.2765
    2.2765   10.8952];
Px22_MUB=[ 4.5393    3.7835
    3.7835    4.7533];


[xx,yy]=meshgrid(-20:0.1:20);
pp1MUB=zeros(size(xx));
pp2MUB=zeros(size(xx));
pp1MI4=zeros(size(xx));
pp2MI4=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        pp1MUB(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px11_MUB);
        pp2MUB(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px22_MUB);
        pp1MI4(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px11_MI4);
        pp2MI4(i,j)=mvnpdf([xx(i,j),yy(i,j)],[0,0],Px22_MI4);
     end
end
figure(1)
contour(xx,yy,pp1MUB,1e-3,'r','linewidth',2)
hold on
contour(xx,yy,pp1MI4,1e-3,'b.-','linewidth',2)
legend('MI_{UB}','MI_4')
xlabel('x')
ylabel('y')
plot_prop_paper
grid


figure(2)
contour(xx,yy,pp2MUB,1e-3,'r','linewidth',2)
hold on
contour(xx,yy,pp2MI4,1e-3,'b.-','linewidth',2)
legend('MI_{UB}','MI_4')
xlabel('x')
ylabel('y')
plot_prop_paper
grid

