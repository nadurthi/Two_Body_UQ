
  matlabpool open 2
 r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
%r2=[-6759.23773137801 ,2109.13744737105,287.210078849733,  0.915903262553281,0.720346730782246 , -7.3978787796829];
% r2=[6729.43097302094	-318.159560024553	-1381.89229666774	2.26490482713380	-1.18931778372278	7.16940625613625]';                  
%T2=0:30:21600;
r2=[ -4849.50137253495,-244.375066964743,5390.17987394301,-4.82791934267752,1.25847522458469 ,-5.1905381510445];

opt = odeset('reltol',1e-12,'abstol',1e-12);

N=300;
% T1=0:30:160740;
% T2=0:30:160740;
T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,160690,10),[160691:160790]];

P1=blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
P2 =blkdiag(0.001,0.001,0.001,1e-7,1e-7,1e-7);
X2mc01=mvnrnd(r1,P1,N);
X2mc02=mvnrnd(r2,P2,N);

X1mc=zeros(length(T1),6,N);
X2mc=zeros(length(T2),6,N);
parfor i=1:N
    i
   [t,X1]= ode45(@twoBody,T1,X2mc01(i,:)',opt); 
   [t,X2]= ode45(@twoBody,T2,X2mc02(i,:)',opt);
   X1mc(:,:,i)=X1;
   X2mc(:,:,i)=X2;
end
tic
[t,X1]= ode45(@twoBody,linspace(0,160791,30),X2mc01(i,:)',opt); 
toc
N
    
  % all done, close pool of workers
  matlabpool close

  
  
  for tk=40:1:110
  XX1=zeros(N,6);
  XX2=zeros(N,6);
  for i=1:1:N
      XX1(i,:)=X1mc(tk,:,i);
      XX2(i,:)=X2mc(tk,:,i);
  end
%   plot3(XX1(:,1),XX1(:,2),XX1(:,3),'bo',XX2(:,1),XX2(:,2),XX2(:,3),'ro')
subplot(1,3,1)

plot(XX1(:,1),XX1(:,2),'bo',XX2(:,1),XX2(:,2),'ro')
subplot(1,3,2)
plot(XX1(:,2),XX1(:,3),'bo',XX2(:,2),XX2(:,3),'ro')
subplot(1,3,3)
plot(XX1(:,1),XX1(:,3),'bo',XX2(:,1),XX2(:,3),'ro')
title(num2str(tk))
pause
  end
  
  
  