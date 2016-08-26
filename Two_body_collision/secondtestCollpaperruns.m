
r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
r2=[6729.43097302094	-318.159560024553	-1381.89229666774	2.26490482713380	-1.18931778372278	7.16940625613625]';


T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];

P1=blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
P2 =blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);

%%

    Nmc=1;
    X1mc0=mvnrnd(r1,P1,Nmc);
    X2mc0=mvnrnd(r2,P2,Nmc);
    
    X1mc=zeros(length(T1),6,Nmc);
    X2mc=zeros(length(T2),6,Nmc);
    
    parfor i=1:Nmc
        i
        tic
        [t,X1]= ode45(@twoBody,T1,X1mc0(i,:)',opt);
        toc
        [t,X2]= ode45(@twoBody,T2,X2mc0(i,:)',opt);
        X1mc(:,:,i)=X1;
        X2mc(:,:,i)=X2;
    end
    w1mc=ones(Nmc,1)/Nmc;
    w2mc=ones(Nmc,1)/Nmc;

	c=clock;
	c(end)=round(c(end)*1000);
	c=num2str(c);
	c(c==' ')='';
	strg=strcat('CollpaperSimsMC_',c);
    save(strg,'X1mc','X2mc','w1mc','w2mc','r1','r2','P1','P2','T1','T2','-v7.3')
    for j=1:1:1000
        if exist(strcat('CollpaperSimsMC_',num2str(j),'.mat'),'file')==0
            break
        end
    end
    %movefile(strcat(strg,'.mat'),strcat('CollpaperSimsMC_',num2str(j),'.mat'))

    %save(strcat('CollpaperSimsMC_',num2str(j)),'X1mc','X2mc','w1mc','w2mc','r1','r2','P1','P2','T1','T2','-v7.3')




%%
160790-21650
  [t,X2]= ode45(@twoBody,T2,X2mc0(i,:)',opt);