load('body2simsMC_final')
load('body2simsOTHERS_final')

time.t0 = 0;
time.dt = 29.173785213265145;
time.tf = 21600;
T=time.t0:time.dt:time.tf;

t=T/3600;
figure
plot(t,sqrt(sum((Mmc1_1e4-Mmc1).^2,2)),'k--',t,sqrt(sum((Mut1-Mmc1).^2,2)),t,sqrt(sum((Mcut41-Mmc1).^2,2)),t,sqrt(sum((Mcut61-Mmc1).^2,2)),t,sqrt(sum((Mcut81-Mmc1).^2,2)),t,sqrt(sum((Mgh31-Mmc1).^2,2)),t,sqrt(sum((Mgh41-Mmc1).^2,2)),t,sqrt(sum((Mgh51-Mmc1).^2,2)),'linewidth',2)
legend('MC-10000','UT','CUT4','CUT6','CUT8','GH3','GH4','GH5')

figure
plot(t,sqrt(sum((Mmc2_1e4-Mmc2).^2,2)),'k--',t,sqrt(sum((Mut2-Mmc2).^2,2)),t,sqrt(sum((Mcut42-Mmc2).^2,2)),t,sqrt(sum((Mcut62-Mmc2).^2,2)),t,sqrt(sum((Mcut82-Mmc2).^2,2)),t,sqrt(sum((Mgh32-Mmc2).^2,2)),t,sqrt(sum((Mgh42-Mmc2).^2,2)),t,sqrt(sum((Mgh52-Mmc2).^2,2)),'linewidth',2)
legend('MC-10000','UT','CUT4','CUT6','CUT8','GH3','GH4','GH5')

figure
plot(t,sqrt(sum((Mmc3_1e4-Mmc3).^2,2)),'k--',t,sqrt(sum((Mut3-Mmc3).^2,2)),t,sqrt(sum((Mcut43-Mmc3).^2,2)),t,sqrt(sum((Mcut63-Mmc3).^2,2)),t,sqrt(sum((Mcut83-Mmc3).^2,2)),t,sqrt(sum((Mgh33-Mmc3).^2,2)),t,sqrt(sum((Mgh43-Mmc3).^2,2)),t,sqrt(sum((Mgh53-Mmc3).^2,2)),'linewidth',2)
legend('MC-10000','UT','CUT4','CUT6','CUT8','GH3','GH4','GH5')

figure
plot(t,sqrt(sum((Mmc4_1e4-Mmc4).^2,2)),'k--',t,sqrt(sum((Mut4-Mmc4).^2,2)),t,sqrt(sum((Mcut44-Mmc4).^2,2)),t,sqrt(sum((Mcut64-Mmc4).^2,2)),t,sqrt(sum((Mcut84-Mmc4).^2,2)),t,sqrt(sum((Mgh34-Mmc4).^2,2)),t,sqrt(sum((Mgh44-Mmc4).^2,2)),t,sqrt(sum((Mgh54-Mmc4).^2,2)),'linewidth',2)
legend('MC-10000','UT','CUT4','CUT6','CUT8','GH3','GH4','GH5')
