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
endtime.min=50;
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
 
N=1000;

P0 =blkdiag(1e-11,1e-11,1e-11,1e-5,1e-5,1e-5);
XmcIR0=mvnrnd([satrecIR.a,satrecIR.ecco,satrecIR.inclo,satrecIR.argpo,satrecIR.nodeo,satrecIR.mo],P0,N);
XmcCM0=mvnrnd([satrecCM.a,satrecCM.ecco,satrecCM.inclo,satrecCM.argpo,satrecCM.nodeo,satrecCM.mo],P0,N);
w=1/N*ones(N,1);

T=0:1:500;
% XIR=zeros(length(T),6);
% XCM=zeros(length(T),6);
% d=zeros(length(T),1);
% for t =1:1:length(T)
%     
%     [satrecIRmc, rIR, vIR] = sgp4(satrecIR,stopmfeIR+T(t));
%     [satrecCMmc, rCM, vCM] = sgp4(satrecCM,stopmfeCM+T(t));
%     XIR(t,:)=[rIR(:)',vIR(:)'];
%     XCM(t,:)=[rCM(:)',vCM(:)'];
%     d(t)=norm(rIR(:)-rCM(:));
% end




XmcIR=zeros(length(T),6,N);
XmcCM=zeros(length(T),6,N);

for i=1:1:N
 i
    satrecIRmc =sim_satrec_orbspace(satrecIR,XmcIR0(i,1),XmcIR0(i,2),XmcIR0(i,3),XmcIR0(i,4),XmcIR0(i,5),XmcIR0(i,6));
   satrecCMmc =sim_satrec_orbspace(satrecCM,XmcCM0(i,1),XmcCM0(i,2),XmcCM0(i,3),XmcCM0(i,4),XmcCM0(i,5),XmcCM0(i,6));
for t =1:1:length(T)
    
    [satrecIRmc, rIR, vIR] = sgp4(satrecIRmc,stopmfeIR+T(t));
    [satrecCMmc, rCM, vCM] = sgp4(satrecCMmc,stopmfeCM+T(t));
    XmcIR(t,:,i)=[rIR(:)',vIR(:)'];
    XmcCM(t,:,i)=[rCM(:)',vCM(:)'];
end
end
keyboard

ptxyz='C:\Users\Nagnanamus\Google Drive\2 body problem\xyzfigs\';
ptxy='C:\Users\Nagnanamus\Google Drive\2 body problem\xyfigs\';
ptyz='C:\Users\Nagnanamus\Google Drive\2 body problem\yzfigs\';
ptxz='C:\Users\Nagnanamus\Google Drive\2 body problem\xzfigs\';

nFrames = length(T);
% Preallocate movie structure.
movxyz(1:nFrames) = struct('cdata', [],'colormap', []);
movxy(1:nFrames) = struct('cdata', [],'colormap', []);
movyz(1:nFrames) = struct('cdata', [],'colormap', []);
movxz(1:nFrames) = struct('cdata', [],'colormap', []);
[xx,yy,zz] = sphere;
for t=1:1:length(T)
XIR=zeros(N,3);
XCM=zeros(N,3);
for i=1:1:N
    XIR(i,:)=XmcIR(t,1:3,i);
    XCM(i,:)=XmcCM(t,1:3,i);
end
[year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg((stopmfeIR+T(t))/1440+satrecIR.jdsatepoch) ;
figure(1)
% subplot(2,2,1)
plot3(XIR(:,1),XIR(:,2),XIR(:,3),'ro',XCM(:,1),XCM(:,2),XCM(:,3),'bo')
hold on
surf(6380.72609981813*xx,6380.72609981813*yy,6380.72609981813*zz,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
colormap(rgb2hsv([0.3 0.9 0.99]))
alpha 0.25
xlabel('x')
ylabel('y')
zlabel('z')
legend('IR','CM')
grid on
title(strcat(num2str(month),'/',num2str(day),'/',num2str(year),'  at ',num2str(hour),':',num2str(minu),':',num2str(sec)),'Interpreter','none')
axis([-8000,8000,-8000,8000,-8000,8000])
% pause(0.5)
movxyz(t) = getframe(gcf);
pause(0.1)
plot_prop_paper
saveas(gcf,strcat(ptxyz,'xyz__',num2str(month),'-',num2str(day),'-',num2str(year),'__',num2str(hour),'-',num2str(minu),'-',num2str(sec)), 'png')
pause(0.1)
close

figure(2)
% subplot(2,2,2)
plot(XIR(:,1),XIR(:,2),'ro',XCM(:,1),XCM(:,2),'bo')
legend('IR','CM')
xlabel('x')
ylabel('y')
grid on
% axis([-5999,1500,-3000,2500])
% axis([-8000,8000,-8000,8000])
title(strcat(num2str(month),'/',num2str(day),'/',num2str(year),' at ',num2str(hour),':',num2str(minu),':',num2str(sec)),'Interpreter','none')
% pause(0.5) 
movxy(t) = getframe(gcf);
pause(0.1)
plot_prop_paper
saveas(gcf,strcat(ptxy,'xy__',num2str(month),'-',num2str(day),'-',num2str(year),'__',num2str(hour),'-',num2str(minu),'-',num2str(sec)), 'png')
pause(0.1)
close

figure(3)
% subplot(2,2,3)
plot(XIR(:,2),XIR(:,3),'ro',XCM(:,2),XCM(:,3),'bo')
legend('IR','CM')
xlabel('y')
ylabel('z')
grid on
% axis([-3000,2500,5000,7500])
title(strcat(num2str(month),'/',num2str(day),'/',num2str(year),' at ',num2str(hour),':',num2str(minu),':',num2str(sec)),'Interpreter','none')
% axis([-8000,8000,-8000,8000])
% pause(0.5)
movyz(t) = getframe(gcf);
pause(0.1)
plot_prop_paper
saveas(gcf,strcat(ptyz,'yz__',num2str(month),'-',num2str(day),'-',num2str(year),'__',num2str(hour),'-',num2str(minu),'-',num2str(sec)), 'png')
pause(0.1)
close

figure(4)
% subplot(2,2,4)
plot(XIR(:,1),XIR(:,3),'ro',XCM(:,1),XCM(:,3),'bo')
legend('IR','CM')
xlabel('x')
ylabel('z')
title(strcat(num2str(month),'/',num2str(day),'/',num2str(year),' at ',num2str(hour),':',num2str(minu),':',num2str(sec)),'Interpreter','none')
% axis([-6000,1500,5000,7500])
% axis([-8000,8000,-8000,8000])
grid on
% pause(0.5)
movxz(t) = getframe(gcf);
pause(0.1)
plot_prop_paper
saveas(gcf,strcat(ptxz,'xz__',num2str(month),'-',num2str(day),'-',num2str(year),'__',num2str(hour),'-',num2str(minu),'-',num2str(sec)), 'png')
pause(0.1)
close


end
 movie(movxyz)
 movie2avi(movxyz, 'MOVxyz.avi', 'compression', 'None');
% movie2avi(movxy, 'MOVxy.avi', 'compression', 'None');
% movie2avi(movyz, 'MOVyz.avi', 'compression', 'None');
% movie2avi(movxz, 'MOVxz.avi', 'compression', 'None');

for t=1:1:length(T)
[year,month,day,hour,minu,sec,dayweek,dategreg] = julian2greg((stopmfeIR+T(t))/1440+satrecIR.jdsatepoch) ;
t
disp(strcat(num2str(month),'/',num2str(day),'/',num2str(year),'  at ',num2str(hour),':',num2str(minu),':',num2str(sec)))

end
XIR=zeros(N,3);
XCM=zeros(N,3);
tt=350;
for i=1:1:N
    XIR(i,:)=XmcIR(tt,1:3,i);
    XCM(i,:)=XmcCM(tt,1:3,i);
end
plot3(XIR(:,1),XIR(:,2),XIR(:,3),'bo',XCM(:,1),XCM(:,2),XCM(:,3),'ro')


plot(XIR(:,2),XIR(:,3),'ro')
plot(XCM(:,2),XCM(:,3),'bo')

plot3(XCM(:,1),XCM(:,2),XCM(:,3),'bo')
plot3(XIR(:,1),XIR(:,2),XIR(:,3),'ro')
grid