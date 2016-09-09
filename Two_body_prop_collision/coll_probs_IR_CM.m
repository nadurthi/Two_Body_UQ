%testing collision probabilities for iridium and cosmos
load('IR-CM_collision_runs')
D=zeros(1001,1);
for tt=1:1:501;
N=size(XmcIR,3);
% N=5000;
XIR=zeros(N,3);
XCM=zeros(N,3);
for i=1:1:N
    XIR(i,:)=XmcIR(tt,1:3,i);
    XCM(i,:)=XmcCM(tt,1:3,i);
end
% plot3(XIR(:,1),XIR(:,2),XIR(:,3),'bo',XCM(:,1),XCM(:,2),XCM(:,3),'ro')
figure(1)
plot(XIR(:,1),XIR(:,2),'bo',XCM(:,1),XCM(:,2),'ro')
% saveas(gca,strcat('xy_',num2str(tt)),'png')
figure(2)
plot(XIR(:,2),XIR(:,3),'bo',XCM(:,2),XCM(:,3),'ro')
% saveas(gca,strcat('yz_',num2str(tt)),'png')
figure(3)
plot(XIR(:,1),XIR(:,3),'bo',XCM(:,1),XCM(:,3),'ro')
% saveas(gca,strcat('xz_',num2str(tt)),'png')
pause(0.01)
%% calculating the 1-1 MC miss distance
d=zeros(N,1);
ind=randperm(N);
k=1;
for i=1:1:N
  
    d(i)=norm(XIR(ind(i),1:3)-XCM(i,1:3));

end
D(tt)=min(d);
% % hist(d)
% pd = fitdist(d,'Gamma');
% x=0:0.01:5;
% pf = pdf(pd,x);
% figure
% plot(x,pf)
end

%% iridium cosmos collision conjuction analyisis
%% From the following TLE data (latest before the collision)
% 2/09/2009 at 11:57:36.8904960001419
CML1='1 22675U 93036A   09040.49834364 -.00000001  00000-0  95251-5 0  7411';
CML2='2 22675  74.0355  19.4646 0016027  98.7014 261.5952 14.31135643817415     0.00      4320.0        360.00';

% 2/09/2009 at 18:49:39.2819519997465  
IRL1='1 24946U 97051C   09040.78448243  .00000153  00000-0  47668-4 0  4775';
IRL2='2 24946  86.3994 121.7028 0002288  85.1644 274.9812 14.34219863597336     0.00      4320.0        360.00' ;

global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
   global opsmode
   opsmode= 'a';
global idebug dbgfile

endtime.year=2009;
endtime.month=2;
endtime.day=10;
endtime.hr=16;
endtime.min=01;
endtime.sec=01;
dtmin=0.5;

[satrecIR, startmfeIR, stopmfeIR, dtmin] = twoline2rv_modified(72,IRL1,IRL2,'m','e',endtime,dtmin);
[satrecCM, startmfeCM, stopmfeCM, dtmin] = twoline2rv_modified(72,CML1,CML2,'m','e',endtime,dtmin);
 [satrecIR, rIR0, vIR0] = sgp4(satrecIR,0);
 [satrecCM, rCM0, vCM0] = sgp4(satrecCM,0);
  
%  [aIR,eIR,iIR,omgIR,OmgIR,MIR]=cart2orb(rIR0(1),rIR0(2),rIR0(3),vIR0(1),vIR0(2),vIR0(3));
%  satrecIRmc = sim_satrec(satrecIR,aIR,eIR,iIR,omgIR,OmgIR,MIR);
%  [satrecIRmc, rIR, vIR] = sgp4(satrecIRmc,0);
%  [rIR0;rIR]
%   [vIR0;vIR]
 
% satrecIRmc =sim_satrec_orbspace(satrecIR,XmcIR0(i,1),XmcIR0(i,2),XmcIR0(i,3),XmcIR0(i,4),XmcIR0(i,5),XmcIR0(i,6));
% satrecCMmc =sim_satrec_orbspace(satrecCM,XmcCM0(i,1),XmcCM0(i,2),XmcCM0(i,3),XmcCM0(i,4),XmcCM0(i,5),XmcCM0(i,6));

T=0:0.167:60;
XIR=zeros(length(T),6);
XCM=zeros(length(T),6);


for t =1:1:length(T)
    
    [satrecIRmc, rIR, vIR] = sgp4(satrecIRmc,stopmfeIR+T(t));
    [satrecCMmc, rCM, vCM] = sgp4(satrecCMmc,stopmfeCM+T(t));
    XIR(t,:)=[rIR(:)',vIR(:)'];
    XCM(t,:)=[rCM(:)',vCM(:)'];
    plot3(XIR(1:t,1),XIR(1:t,2),XIR(1:t,3),'ro',XCM(1:t,1),XCM(1:t,2),XCM(1:t,3),'bo')
pause(0.02)
end
