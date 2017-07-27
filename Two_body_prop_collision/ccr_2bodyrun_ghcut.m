function ccr_2bodyrun_ghcut()
r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
%r2=[-6759.23773137801 ,2109.13744737105,287.210078849733,  0.915903262553281,0.720346730782246 , -7.3978787796829];
r2=[6729.43097302094	-318.159560024553	-1381.89229666774	2.26490482713380	-1.18931778372278	7.16940625613625];                  

matlabpool open 12

opt = odeset('reltol',1e-12,'abstol',1e-12);

% N=100000;
%T1=0:30:160740;
%T2=0:30:21600;
T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];


P1=blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
P2 =blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);
[X01gh4,w1gh4]=GH_points(r1,P1,4);
[X02gh4,w2gh4]=GH_points(r2,P2,4);
N=length(w1gh4);
X1gh4=zeros(length(T1),6,N);
X2gh4=zeros(length(T2),6,N);
parfor i=1:N
    i
   [t,X1]= ode45(@twoBody,T1,X01gh4(i,:)',opt); 
   [t,X2]= ode45(@twoBody,T2,X02gh4(i,:)',opt);
   X1gh4(:,:,i)=X1;
   X2gh4(:,:,i)=X2;
end

[X01cut6,w1cut6]=conjugate_dir_gausspts_till_6moment_scheme2(r1,P1);
[X02cut6,w2cut6]=conjugate_dir_gausspts_till_6moment_scheme2(r2,P2);
N=length(w1cut6);
X1cut6=zeros(length(T1),6,N);
X2cut6=zeros(length(T2),6,N);
parfor i=1:N
    i
   [t,X1]= ode45(@twoBody,T1,X01cut6(i,:)',opt); 
   [t,X2]= ode45(@twoBody,T2,X02cut6(i,:)',opt);
   X1cut6(:,:,i)=X1;
   X2cut6(:,:,i)=X2;
end

[X01cut8,w1cut8]=conjugate_dir_gausspts_till_8moment(r1,P1);
[X02cut8,w2cut8]=conjugate_dir_gausspts_till_8moment(r2,P2);
N=length(w1cut8);
X1cut8=zeros(length(T1),6,N);
X2cut8=zeros(length(T2),6,N);
parfor i=1:N
    i
   [t,X1]= ode45(@twoBody,T1,X01cut8(i,:)',opt); 
   [t,X2]= ode45(@twoBody,T2,X02cut8(i,:)',opt);
   X1cut8(:,:,i)=X1;
   X2cut8(:,:,i)=X2;
end
save('Collisioneg2_GH4_CUT6_CUT8','X1cut8','w1cut8','X2cut8','w2cut8','X1cut6','w1cut6','X2cut6','w2cut6','X1gh5','w1gh5','X2gh5','w2gh5','T','-v7.3')

matlabpool close


end