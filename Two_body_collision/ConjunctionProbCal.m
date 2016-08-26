%% momnet convergence


[y1,M11gh4]=Evol_moments_samples(X1gh4,w1gh4,1,'central');
[y1,M21gh4]=Evol_moments_samples(X2gh4,w2gh4,1,'central');

[y2,M12gh4]=Evol_moments_samples(X1gh4,w1gh4,2,'central');
[y2,M22gh4]=Evol_moments_samples(X2gh4,w2gh4,2,'central');

[y3,M13gh4]=Evol_moments_samples(X1gh4,w1gh4,3,'central');
[y3,M23gh4]=Evol_moments_samples(X2gh4,w2gh4,3,'central');

[y4,M14gh4]=Evol_moments_samples(X1gh4,w1gh4,4,'central');
[y4,M24gh4]=Evol_moments_samples(X2gh4,w2gh4,4,'central');

% %
[y1,M11cut6]=Evol_moments_samples(X1cut6,w1cut6,1,'central');
[y1,M21cut6]=Evol_moments_samples(X2cut6,w2cut6,1,'central');

[y2,M12cut6]=Evol_moments_samples(X1cut6,w1cut6,2,'central');
[y2,M22cut6]=Evol_moments_samples(X2cut6,w2cut6,2,'central');

[y3,M13cut6]=Evol_moments_samples(X1cut6,w1cut6,3,'central');
[y3,M23cut6]=Evol_moments_samples(X2cut6,w2cut6,3,'central');

[y4,M14cut6]=Evol_moments_samples(X1cut6,w1cut6,4,'central');
[y4,M24cut6]=Evol_moments_samples(X2cut6,w2cut6,4,'central');

% %
[y1,M11cut8]=Evol_moments_samples(X1cut8,w1cut8,1,'central');
[y1,M21cut8]=Evol_moments_samples(X2cut8,w2cut8,1,'central');

[y2,M12cut8]=Evol_moments_samples(X1cut8,w1cut8,2,'central');
[y2,M22cut8]=Evol_moments_samples(X2cut8,w2cut8,2,'central');

[y3,M13cut8]=Evol_moments_samples(X1cut8,w1cut8,3,'central');
[y3,M23cut8]=Evol_moments_samples(X2cut8,w2cut8,3,'central');

[y4,M14cut8]=Evol_moments_samples(X1cut8,w1cut8,4,'central');
[y4,M24cut8]=Evol_moments_samples(X2cut8,w2cut8,4,'central');

% %
w1mc=ones(1e5,1)/1e5;
w2mc=w1mc;
[y1,M11mc]=Evol_moments_samples(X1mc,w1mc,1,'central');
[y1,M21mc]=Evol_moments_samples(X2mc,w2mc,1,'central');

[y2,M12mc]=Evol_moments_samples(X1mc,w1mc,2,'central');
[y2,M22mc]=Evol_moments_samples(X2mc,w2mc,2,'central');

[y3,M13mc]=Evol_moments_samples(X1mc,w1mc,3,'central');
[y3,M23mc]=Evol_moments_samples(X2mc,w2mc,3,'central');

[y4,M14mc]=Evol_moments_samples(X1mc,w1mc,4,'central');
[y4,M24mc]=Evol_moments_samples(X2mc,w2mc,4,'central');

save('CollisionBodiesGH_allMOMS','M11gh4','M21gh4','M12gh4','M22gh4','M13gh4','M23gh4','M14gh4','M24gh4','M11cut6','M21cut6','M12cut6','M22cut6','M13cut6','M23cut6','M14cut6','M24cut6','M11cut8','M21cut8','M12cut8','M22cut8','M13cut8','M23cut8','M14cut8','M24cut8','M11mc','M21mc','M12mc','M22mc','M13mc','M23mc','M14mc','M24mc','-append')

%% Convergence of momnets

e11gh4=errovertime(M11gh4,M11gh7)
e21gh4=errovertime(M21gh4,M21gh7)

e12gh4=errovertime(M12gh4,M12gh7)
e22gh4=errovertime(M22gh4,M22gh7)

e13gh4=errovertime(M13gh4,M13gh7)
e23gh4=errovertime(M23gh4,M23gh7)

e14gh4=errovertime(M14gh4,M14gh7)
e24gh4=errovertime(M24gh4,M24gh7)

% %
e11gh5=errovertime(M11gh5,M11gh7)
e21gh5=errovertime(M21gh5,M21gh7)

e12gh5=errovertime(M12gh5,M12gh7)
e22gh5=errovertime(M22gh5,M22gh7)

e13gh5=errovertime(M13gh5,M13gh7)
e23gh5=errovertime(M23gh5,M23gh7)

e14gh5=errovertime(M14gh5,M14gh7)
e24gh5=errovertime(M24gh5,M24gh7)

% %
e11gh6=errovertime(M11gh6,M11gh7)
e21gh6=errovertime(M21gh6,M21gh7)

e12gh6=errovertime(M12gh6,M12gh7)
e22gh6=errovertime(M22gh6,M22gh7)

e13gh6=errovertime(M13gh6,M13gh7)
e23gh6=errovertime(M23gh6,M23gh7)

e14gh6=errovertime(M14gh6,M14gh7)
e24gh6=errovertime(M24gh6,M24gh7)

% %
e11mc=errovertime(M11mc,M11gh7)
e21mc=errovertime(M21mc,M21gh7)

e12mc=errovertime(M12mc,M12gh7)
e22mc=errovertime(M22mc,M22gh7)

e13mc=errovertime(M13mc,M13gh7)
e23mc=errovertime(M23mc,M23gh7)

e14mc=errovertime(M14mc,M14gh7)
e24mc=errovertime(M24mc,M24gh7)

% %
e11cut6=errovertime(M11cut6,M11gh7)
e21cut6=errovertime(M21cut6,M21gh7)

e12cut6=errovertime(M12cut6,M12gh7)
e22cut6=errovertime(M22cut6,M22gh7)

e13cut6=errovertime(M13cut6,M13gh7)
e23cut6=errovertime(M23cut6,M23gh7)

e14cut6=errovertime(M14cut6,M14gh7)
e24cut6=errovertime(M24cut6,M24gh7)

% %
e11cut8=errovertime(M11cut8,M11gh7)
e21cut8=errovertime(M21cut8,M21gh7)

e12cut8=errovertime(M12cut8,M12gh7)
e22cut8=errovertime(M22cut8,M22gh7)

e13cut8=errovertime(M13cut8,M13gh7)
e23cut8=errovertime(M23cut8,M23gh7)

e14cut8=errovertime(M14cut8,M14gh7)
e24cut8=errovertime(M24cut8,M24gh7)


[e11gh4,e11gh5,e11gh6,e11cut6,e11cut8,e11mc;...
 e12gh4,e12gh5,e12gh6,e12cut6,e12cut8,e12mc;...
 e13gh4,e13gh5,e13gh6,e13cut6,e13cut8,e13mc;...
 e14gh4,e14gh5,e14gh6,e14cut6,e14cut8,e14mc;]

[e21gh4,e21gh5,e21gh6,e21cut6,e21cut8,e21mc;...
 e22gh4,e22gh5,e22gh6,e22cut6,e22cut8,e22mc;...
 e23gh4,e23gh5,e23gh6,e23cut6,e23cut8,e23mc;...
 e24gh4,e24gh5,e24gh6,e24cut6,e24cut8,e24mc;]

plot([4^6,5^6,6^6],[e11gh4,e11gh5,e11gh6],'o-','linewidth',2)
plot([4^6,5^6,6^6],[e12gh4,e12gh5,e12gh6],'o-','linewidth',2)
plot([4^6,5^6,6^6],[e13gh4,e13gh5,e13gh6],'o-','linewidth',2)
plot([4^6,5^6,6^6],[e14gh4,e14gh5,e14gh6],'o-','linewidth',2)

plot([4^6,5^6,6^6],[e24gh4,e24gh5,e24gh6],'o-','linewidth',2)
%% miss ditance and plots of moment convergence in miss distance
N=size(X1mc,3);
Nt=size(X1mc,1);
X11mc=zeros(N,3);
X22mc=zeros(N,3);
X11cut8=zeros(size(X1cut8,3),3);
X22cut8=zeros(size(X1cut8,3),3);
X11cut6=zeros(size(X1cut6,3),3);
X22cut6=zeros(size(X1cut6,3),3);
X11gh4=zeros(size(X1gh4,3),3);
X22gh4=zeros(size(X1gh4,3),3);
X11gh5=zeros(size(X1gh5,3),3);
X22gh5=zeros(size(X1gh5,3),3);
X11gh6=zeros(size(X1gh6,3),3);
X22gh6=zeros(size(X1gh6,3),3);
X11gh7=zeros(size(X1gh7,3),3);
X22gh7=zeros(size(X1gh7,3),3);
t=40:70;
D=cell(length(t),1);
Dmc=cell(length(t),1);
Dcut8=cell(length(t),1);
Dcut6=cell(length(t),1);
Dgh4=cell(length(t),1);
Dgh5=cell(length(t),1);
Dgh6=cell(length(t),1);
Dgh7=cell(length(t),1);
k=1;
for i=t 
%     D{k}=zeros(N,3);
%     tic
%     for j=1:1:N
%         X11mc(j,:)=X1mc(i,1:3,j);
%         X22mc(j,:)=X2mc(i,1:3,j);
%         D{k}(j,:)=X1(j,:)-X2(j,:);
%     end
%     Dmc{k}=missdistpts(X11mc,w1mc,X22mc,w2mc);
%     toc
% tic
%     for j=1:1:size(X1cut8,3)
%         X11cut8(j,:)=X1cut8(i,1:3,j);
%         X22cut8(j,:)=X2cut8(i,1:3,j);
%     end
%     Dcut8{k}=missdistpts(X11cut8,w1cut8,X22cut8,w2cut8);
% toc    
%     for j=1:1:size(X1cut6,3)
%         X11cut6(j,:)=X1cut6(i,1:3,j);
%         X22cut6(j,:)=X2cut6(i,1:3,j);
%     end
%     Dcut6{k}=missdistpts(X11cut6,w1cut6,X22cut6,w2cut6);
% tic    
%     for j=1:1:size(X1gh4,3)
%         X11gh4(j,:)=X1gh4(i,1:3,j);
%         X22gh4(j,:)=X2gh4(i,1:3,j);
%     end
%         Dgh4{k}=missdistpts(X11gh4,w1gh4,X22gh4,w2gh4);
%  toc 
 tic 
    for j=1:1:size(X1gh5,3)
        X11gh5(j,:)=X1gh5(i,1:3,j);
        X22gh5(j,:)=X2gh5(i,1:3,j);
    end
        Dgh5{k}=missdistpts(X11gh5,w1gh5,X22gh5,w2gh5);
 toc
%  tic 
%     for j=1:1:size(X1gh6,3)
%         X11gh6(j,:)=X1gh6(i,1:3,j);
%         X22gh6(j,:)=X2gh6(i,1:3,j);
%     end
%     Dgh6{k}=missdistpts(X11gh6,w1gh6,X22gh6,w2gh6);
% toc
%     for j=1:1:size(X1gh7,3)
%         X11gh7(j,:)=X1gh7(i,1:3,j);
%         X22gh7(j,:)=X2gh7(i,1:3,j);
%     end
%     Dgh7{k}=missdistpts(X11gh7,w1gh7,X22gh7,w2gh7);
%     
    k=k+1
end    
save('MissDistanceMoms','Dgh5','-append')

MisserrM1gh4=0;
MisserrM1gh5=0;
MisserrM1cut6=0;
MisserrM1cut8=0;

MisserrM2gh4=0;
MisserrM2gh5=0;
MisserrM2cut6=0;
MisserrM2cut8=0;

MisserrM3gh4=0;
MisserrM3gh5=0;
MisserrM3cut6=0;
MisserrM3cut8=0;

MisserrM4gh4=0;
MisserrM4gh5=0;
MisserrM4cut6=0;
MisserrM4cut8=0;

for k=1:1:16
MisserrM1gh4(k)=sqrt(sum((Dgh4{k}.M1-Dgh6{k}.M1).^2)/3);
MisserrM1gh5(k)=sqrt(sum((Dgh5{k}.M1-Dgh6{k}.M1).^2)/3);
MisserrM1cut6(k)=sqrt(sum((Dcut6{k}.M1-Dgh6{k}.M1).^2)/3);
MisserrM1cut8(k)=sqrt(sum((Dcut8{k}.M1-Dgh6{k}.M1).^2)/3);

MisserrM2gh4(k)=sqrt(sum((Dgh4{k}.M2-Dgh6{k}.M2).^2)/6);
MisserrM2gh5(k)=sqrt(sum((Dgh5{k}.M2-Dgh6{k}.M2).^2)/6);
MisserrM2cut6(k)=sqrt(sum((Dcut6{k}.M2-Dgh6{k}.M2).^2)/6);
MisserrM2cut8(k)=sqrt(sum((Dcut8{k}.M2-Dgh6{k}.M2).^2)/6);

MisserrM3gh4(k)=sqrt(sum((Dgh4{k}.M3-Dgh6{k}.M3).^2)/10);
MisserrM3gh5(k)=sqrt(sum((Dgh5{k}.M3-Dgh6{k}.M3).^2)/10);
MisserrM3cut6(k)=sqrt(sum((Dcut6{k}.M3-Dgh6{k}.M3).^2)/10);
MisserrM3cut8(k)=sqrt(sum((Dcut8{k}.M3-Dgh6{k}.M3).^2)/10);

MisserrM4gh4(k)=sqrt(sum((Dgh4{k}.M4-Dgh6{k}.M4).^2)/15);
MisserrM4gh5(k)=sqrt(sum((Dgh5{k}.M4-Dgh6{k}.M4).^2)/15);
MisserrM4cut6(k)=sqrt(sum((Dcut6{k}.M4-Dgh6{k}.M4).^2)/15);
MisserrM4cut8(k)=sqrt(sum((Dcut8{k}.M4-Dgh6{k}.M4).^2)/15);

end

[sqrt(sum(MisserrM1gh4.^2)/16),sqrt(sum(MisserrM1gh5.^2)/16),sqrt(sum(MisserrM1cut6.^2)/16),sqrt(sum(MisserrM1cut8.^2)/16)]
[sqrt(sum(MisserrM2gh4.^2)/16),sqrt(sum(MisserrM2gh5.^2)/16),sqrt(sum(MisserrM2cut6.^2)/16),sqrt(sum(MisserrM2cut8.^2)/16)]
[sqrt(sum(MisserrM3gh4.^2)/16),sqrt(sum(MisserrM3gh5.^2)/16),sqrt(sum(MisserrM3cut6.^2)/16),sqrt(sum(MisserrM3cut8.^2)/16)]
[sqrt(sum(MisserrM4gh4.^2)/16),sqrt(sum(MisserrM4gh5.^2)/16),sqrt(sum(MisserrM4cut6.^2)/16),sqrt(sum(MisserrM4cut8.^2)/16)]
%%

load('MissDistanceMoms')


T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];

N=size(X1mc,3);
Nt=size(X1mc,1);
X1=zeros(N,3);
X2=zeros(N,3);
R=[0.01:0.01:0.1,0.2:0.1:0.5];
t=40:70;
probcollcut8=zeros(length(R),length(t));
probcollcut8_5moms=zeros(length(R),length(t));

probcollcut6=zeros(1,length(t));


probcollgh4=zeros(1,length(t));
probcollgh5=zeros(1,length(t));
probcollgh6=zeros(length(R),length(t));
probcollgh6_5moms=zeros(length(R),length(t));

probgauss=zeros(length(R),length(t));

X11cut8=zeros(size(X1cut8,3),6);
X22cut8=zeros(size(X1cut8,3),6);

X11cut6=zeros(size(X1cut6,3),3);
X22cut6=zeros(size(X1cut6,3),3);

X11gh4=zeros(size(X1gh4,3),3);
X22gh4=zeros(size(X1gh4,3),3);

X11gh5=zeros(size(X1gh5,3),3);
X22gh5=zeros(size(X1gh5,3),3);

X11gh6=zeros(size(X1gh6,3),3);
X22gh6=zeros(size(X1gh6,3),3);

k=1;
D=zeros(N,3);
for i=t 
%     for j=1:1:N
%         X1(j,:)=X1mc(i,1:3,j);
%         X2(j,:)=X2mc(i,1:3,j);
%         D(j,:)=X1(j,:)-X2(j,:);
%     end
%     figure(1)
%     plot3(X1(1:1e4,1),X1(1:1e4,2),X1(1:1e4,3),'ro',X2(1:1e4,1),X2(1:1e4,2),X2(1:1e4,3),'bo')
%     view(248,22)
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     grid
%     plot_prop_paper
%     saveas(gcf,strcat('C:\Users\Nagnanamus\Google Drive\Nagavenkat_Adurthi_DRIVE\2014 papers\JGCD space collision 2014\figures\sim',num2str(k)),'pdf')
%     figure(2)
%     plot3(D(1:1e4,1),D(1:1e4,2),D(1:1e4,3),'ko')
%     view(46,38)
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     grid
%     plot_prop_paper
%     saveas(gcf,strcat('C:\Users\Nagnanamus\Google Drive\Nagavenkat_Adurthi_DRIVE\2014 papers\JGCD space collision 2014\figures\missdistsim',num2str(k)),'pdf')
%     pause(1) 
    
    for j=1:1:size(X1cut8,3)
        X11cut8(j,:)=X1cut8(i,1:6,j);
        X22cut8(j,:)=X2cut8(i,1:6,j);
    end
%     for j=1:1:size(X1cut6,3)
%         X11cut6(j,:)=X1cut6(i,1:3,j);
%         X22cut6(j,:)=X2cut6(i,1:3,j);
%     end
%     for j=1:1:size(X1gh4,3)
%         X11gh4(j,:)=X1gh4(i,1:3,j);
%         X22gh4(j,:)=X2gh4(i,1:3,j);
%     end
    for j=1:1:size(X1gh5,3)
        X11gh5(j,:)=X1gh5(i,1:3,j);
        X22gh5(j,:)=X2gh5(i,1:3,j);
    end
    for j=1:1:size(X1gh6,3)
        X11gh6(j,:)=X1gh6(i,1:3,j);
        X22gh6(j,:)=X2gh6(i,1:3,j);
    end
% prob_collmc(k)=samples_prob_coll(X1,X2,R,1);
% probnormgh4(k)=prob_coll_norms(X11gh4,w1gh4,X22gh4,w2gh4,R);
%  probnormcut8(k)=prob_coll_norms(X11cut8,w1cut8,X22cut8,w2cut8,R);
% probnormgh(k)=prob_coll_norms(X11gh5,w1gh5,X22gh5,w2gh5,R);
i  
% plot3(D(:,1),D(:,2),D(:,3),'bo')
% keyboard
  probcollcut8(:,k)=prob_coll_nongauss_modf(Dcut8{i-49},4,R,i);
%   probcollcut8_5moms(:,k)=prob_coll_nongauss_modf(Dcut8{i-49},5,R);
%  probcollcut8_6moms(:,k)=prob_coll_nongauss_modf(Dcut8{i-49},6,R);
 
%  probcollgh6(:,k)=prob_coll_nongauss_modf(Dgh6{i-49},4,R);
%  probcollgh6_5moms(:,k)=prob_coll_nongauss_modf(Dgh6{i-49},5,R);
% probcollcut6(k)=prob_coll_nongauss(X11cut6,w1cut6,X22cut6,w2cut6,R);

% probcollgh4(k)=prob_coll_nongauss(X11gh4,w1gh4,X22gh4,w2gh4,R);
% probcollgh5(k)=prob_coll_nongauss(X11gh5,w1gh5,X22gh5,w2gh5,R);
% probcollgh6(k)=prob_coll_nongauss(X11gh6,w1gh6,X22gh6,w2gh6,R);

%  probnongauss(k)=prob_coll_nongauss(X11gh4,w1gh4,X22gh4,w2gh4,R);
[mu1,P1]=ptswts2muP(X11gh6,w1gh6);
[mu2,P2]=ptswts2muP(X22gh6,w2gh6);
probgauss(:,k)=prob_gaussian_coll(mu1,P1,mu2,P2,R,i);

k=k+1;
end

% save('CUT8GH6probnongausscoll45moms','probgauss','probcollcut8','probcollcut8_5moms','probcollgh6','probcollgh6_5moms')

prob_collmc_square_avg=0;
prob_collmc_circle_avg=0;
for i=1:1:6000
    prob_collmc_square_avg=prob_collmc_square_avg+ProbCollMCSquare6000{i};
    prob_collmc_circle_avg=prob_collmc_circle_avg+ProbCollMCCircle6000{i};
end
prob_collmc_square_avg=prob_collmc_square_avg/6000;
prob_collmc_circle_avg=prob_collmc_circle_avg/6000;

prob_collmc_square_var=0;
prob_collmc_circle_var=0;
for i=1:1:6000
    prob_collmc_square_var=prob_collmc_square_var+(ProbCollMCSquare6000{i}-prob_collmc_square_avg).^2;
    prob_collmc_circle_var=prob_collmc_circle_var+(ProbCollMCCircle6000{i}-prob_collmc_circle_avg).^2;
end
prob_collmc_square_var=prob_collmc_square_var/6000;
prob_collmc_circle_var=prob_collmc_circle_var/6000;

prob_collmc_square_std_ub=prob_collmc_square_avg+1*prob_collmc_square_var.^(0.5);
prob_collmc_square_std_lb=prob_collmc_square_avg-1*prob_collmc_square_var.^(0.5);

for indd=1:1:14
TT=T1(t);
pp=6;
TT=TT(pp:end);
figure(1)
semilogy(TT,prob_collmc_square_avg(indd,pp:end),'o-',TT,probgauss(indd,pp:end),'s-',TT,probcollcut8(indd,pp:end),'^-',TT,probcollcut8_5moms(indd,pp:end),'*-',TT,prob_collmc_square_std_lb(indd,pp:end),'k--',TT,prob_collmc_square_std_ub(indd,pp:end),'k--','linewidth',2,'MarkerSize',12)
set(gca,'XTick',TT)
set(gca,'XTickLabel',{'tca-5','tca-4','tca-3','tca-2','tca-1','tca','tca+1','tca+2','tca+3','tca+4','tca+5'})
legend('MC','Gaussian','CUT8-4','CUT8-5')
title(strcat('R = ',num2str(R(indd)*1000),' (m)'))
xlabel('time in seconds')
ylabel('Probability of collision for threshold R')
plot_prop_paper
saveas(gcf,strcat('CUT8eg1coll',num2str(R(indd)*1000)),'pdf')

figure(2)
semilogy(TT,prob_collmc_square_avg(indd,pp:end),'o-',TT,probgauss(indd,pp:end),'s-',TT,probcollgh6(indd,pp:end),'^-',TT,probcollgh6_5moms(indd,pp:end),'*-',TT,prob_collmc_square_std_lb(indd,pp:end),'k--',TT,prob_collmc_square_std_ub(indd,pp:end),'k--','linewidth',2,'MarkerSize',12)
set(gca,'XTick',TT)
set(gca,'XTickLabel',{'tca-5','tca-4','tca-3','tca-2','tca-1','tca','tca+1','tca+2','tca+3','tca+4','tca+5'})
legend('MC','Gaussian','GH6-4','GH6-5')
title(strcat('R = ',num2str(R(indd)*1000),' (m)'))
xlabel('time in seconds')
ylabel('Probability of collision for threshold R')
plot_prop_paper
saveas(gcf,strcat('GH6eg1coll',num2str(R(indd)*1000)),'pdf')

end
% '+3\sigma_{MC}','-3\sigma_{MC}'
% TT(7:end),prob_collmc_square_std_ub(indd,7:end),'k--',TT(7:end),prob_collmc_square_std_lb(indd,7:end),'k--',

%% 4moms 


load('MissDistanceMoms')
load('CollpaperSimsCUT8')

T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];

R=0.4;
t=40:70;
probcollcut8=zeros(length(t),1);





k=1;

FullRelMotionGaussPdf_mu=cell(length(t),1);
FullRelMotionGaussPdf_P=cell(length(t),1);

filename='CUT8PMEpdfs_4moms';

YY=cell(length(t),1);
parfor i=t 
i
    X11cut8=zeros(size(X1cut8,3),6);
    X22cut8=zeros(size(X1cut8,3),6);
    for j=1:1:size(X1cut8,3)
        X11cut8(j,:)=X1cut8(i,1:6,j);
        X22cut8(j,:)=X2cut8(i,1:6,j);
    end

% keyboard
YY{i-39}=PMEmissPDFS(Dcut8{i-39},4,i-39,X11cut8(:,1:3),X22cut8(:,1:3),w1cut8,w2cut8,'abcd',filename);


[mu1,P1]=ptswts2muP(X11cut8,w1cut8);
[mu2,P2]=ptswts2muP(X22cut8,w2cut8);
mur=mu1-mu2;
Pr=P1+P2;
FullRelMotionGaussPdf_mu{i-39,1}=mur;
FullRelMotionGaussPdf_P{i-39,1}=Pr;

k=k+1;
end


save(filename,'YY','FullRelMotionGaussPdf_mu','FullRelMotionGaussPdf_P')

%% velocity pdfs
load('MissDistanceMoms')
load('CollpaperSimsCUT8')

T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];

R=0.4;
t=40:70;
probcollcut8=zeros(length(t),1);





k=1;

FullRelMotionGaussPdf_mu=cell(length(t),1);
FullRelMotionGaussPdf_P=cell(length(t),1);

filename='CUT8PMEVELpdfs';

YY=cell(length(t),1);
for i=t 
i
    X11cut8=zeros(size(X1cut8,3),6);
    X22cut8=zeros(size(X1cut8,3),6);
    for j=1:1:size(X1cut8,3)
        X11cut8(j,:)=X1cut8(i,1:6,j);
        X22cut8(j,:)=X2cut8(i,1:6,j);
    end

% keyboard
YY{i-39}=PMEmissPDFS(Dcut8{i-39},4,i-39,X11cut8(:,4:6),X22cut8(:,4:6),w1cut8,w2cut8,'abcd',filename);


[mu1,P1]=ptswts2muP(X11cut8,w1cut8);
[mu2,P2]=ptswts2muP(X22cut8,w2cut8);
mur=mu1-mu2;
Pr=P1+P2;
FullRelMotionGaussPdf_mu{i-39,1}=mur;
FullRelMotionGaussPdf_P{i-39,1}=Pr;

k=k+1;
end


save(filename,'YY','FullRelMotionGaussPdf_mu','FullRelMotionGaussPdf_P')

close all
plot_PMEpdfMarg(YY{1},R)


for i=1:1:31
    YY{i}.muu'
    
end
%% Fullpdf


load('MissDistanceMoms')
load('CollpaperSimsCUT8')

T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];

R=0.4;
t=40:70;
probcollcut8=zeros(length(t),1);





k=1;

FullRelMotionGaussPdf_mu=cell(length(t),1);
FullRelMotionGaussPdf_P=cell(length(t),1);

filename='CUT8PMEpdfs';

YY=cell(length(t),1);
for i=t 
i
    X11cut8=zeros(size(X1cut8,3),6);
    X22cut8=zeros(size(X1cut8,3),6);
    for j=1:1:size(X1cut8,3)
        X11cut8(j,:)=X1cut8(i,1:6,j);
        X22cut8(j,:)=X2cut8(i,1:6,j);
    end

% keyboard
YY{i-39}=PMEmissPDFS(Dcut8{i-39},4,i-39,X11cut8(:,1:6),X22cut8(:,1:6),w1cut8,w2cut8,'abcd',filename);


[mu1,P1]=ptswts2muP(X11cut8,w1cut8);
[mu2,P2]=ptswts2muP(X22cut8,w2cut8);
mur=mu1-mu2;
Pr=P1+P2;
FullRelMotionGaussPdf_mu{i-39,1}=mur;
FullRelMotionGaussPdf_P{i-39,1}=Pr;

k=k+1;
end


save(filename,'YY','FullRelMotionGaussPdf_mu','FullRelMotionGaussPdf_P')
%%
load('CUT8PMEpdfs_4moms')
t=40:70;
R=0.4;
F=zeros(length(t),1);
Fg=zeros(length(t),1);

for i=1:1:length(t)
    i
    murv=FullRelMotionGaussPdf_mu{i};
    murv=murv(:);
    PP=FullRelMotionGaussPdf_P{i};
    mur=murv(1:3);
    muv=murv(4:6);
    Pr=PP([1:3],[1:3]);
    Pv=PP([4:6],[4:6]);
    Prv=PP([1:3],[4:6]);
    Ys=cell(1,5);
%     Ys{1}=YY{i,1}; Ys{2}=YY{i,2}; Ys{3}=YY{i,3}; Ys{4}=YY{i,4}; Ys{5}=YY{i,5};
   Ys{1}=YY{i}.Y;
Ys{2}=YY{i}.lam;
Ys{3}=YY{i}.muu;
Ys{4}=YY{i}.iP;
Ys{5}=YY{i}.pdf;

     X11cut8vel=zeros(size(X1cut8,3),3);
    X22cut8vel=zeros(size(X1cut8,3),3);
    for j=1:1:size(X1cut8,3)
        X11cut8vel(j,:)=X1cut8(t(i),4:6,j);
        X22cut8vel(j,:)=X2cut8(t(i),4:6,j);
    end
    
    F(i)=CummCollProb_t(R,mur,muv,Pr,Prv,Pv,Ys,'nongauss',X11cut8vel,X22cut8vel,w1cut8,w2cut8);
    
    % now gaussian
    
    Fg(i)=CummCollProb_t(R,mur,muv,Pr,Prv,Pv,Ys,'gauss',X11cut8vel,X22cut8vel,w1cut8,w2cut8);
end

Q=zeros(length(t),1);
Qg=zeros(length(t),1);

for i=1:1:length(t)
Q(i) = trapz(F(1:i));
Qg(i) = trapz(Fg(1:i));

end

%%
load('CUT8PMEpdfs_5moms')
t=40:70;
R=0.4;
F=zeros(length(t),1);
Fg=zeros(length(t),1);

for i=1:1:length(t)
    i
    murv=FullRelMotionGaussPdf_mu{i};
    murv=murv(:);
    PP=FullRelMotionGaussPdf_P{i};
    mur=murv(1:3);
    muv=murv(4:6);
    Pr=PP([1:3],[1:3]);
    Pv=PP([4:6],[4:6]);
    Prv=PP([1:3],[4:6]);
    Ys=cell(1,5);
    Ys{1}=YY{i,1}; Ys{2}=YY{i,2}; Ys{3}=YY{i,3}; Ys{4}=YY{i,4}; Ys{5}=YY{i,5};
    
    X11cut8vel=zeros(size(X1cut8,3),3);
    X22cut8vel=zeros(size(X1cut8,3),3);
    for j=1:1:size(X1cut8,3)
        X11cut8vel(j,:)=X1cut8(t(i),4:6,j);
        X22cut8vel(j,:)=X2cut8(t(i),4:6,j);
    end
    
    F(i)=CummCollProb_t(R,mur,muv,Pr,Prv,Pv,Ys,'nongauss',X11cut8vel,X22cut8vel,w1cut8,w2cut8);
    
    % now gaussian
    
    Fg(i)=CummCollProb_t(R,mur,muv,Pr,Prv,Pv,Ys,'gauss',X11cut8vel,X22cut8vel,w1cut8,w2cut8);
end

Qmom5=zeros(length(t),1);
Qgmom5=zeros(length(t),1);

for i=1:1:length(t)
Qmom5(i) = trapz(F(1:i));
Qgmom5(i) = trapz(Fg(1:i));

end

%%


load('MCcummProbColl_100')

Pmc=zeros(31,length(ProbcumMCNN));
for i=1:1:length(ProbcumMCNN)
    Pmc(:,i)=ProbcumMCNN{i};
end
Pmc_mu=mean(Pmc,2);
Pmc_std=std(Pmc,0,2);
Pmc_max=max(Pmc,[],2);

figure
plot(t,Q,t,Qg,t,Pmc_mu,t,Pmc_mu+3*Pmc_std,'k--',t,Pmc_mu-3*Pmc_std,'k--')
legend('Non-Gaussian','Gaussian','MC','MC+std','MC-std')

figure
plot(t,Q,t,Qg,t,Pmc_max,t,Pmc_max_20,'linewidth',2)
legend('Non-Gaussian','Gaussian','MC10','MC20')


figure
plot(t,Q,t,Qmom5,t,Qg,t,Pmc_max,t,Pmc_max_20,'linewidth',2)
legend('Non-Gaussian-4','Non-Gaussian-5','Gaussian','MC10','MC20')
set(gca,'YScale','log');

TT=T1(50:70)
plot(TT,Q(11:end),'linewidth',2)
xlabel('time')
ylabel('P_c')
% plot_prop_paper
set(gca,'YTick',linspace(0,2e-4,10))
set(gca,'XTick',TT)
set(gca,'XTickLabel',{'tca-10','tca-9','tca-8','tca-7','tca-6','tca-5','tca-4','tca-3','tca-2','tca-1','tca','tca+1','tca+2','tca+3','tca+4','tca+5','tca+6','tca+7','tca+8','tca+9','tca+10'})
rotateticklabel(gca)
