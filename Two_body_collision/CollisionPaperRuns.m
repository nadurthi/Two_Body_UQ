function CollisionPaperRuns(numWorkers,method)
% start by opening a lab with desired number of workers

% matlabpool('local', numWorkers)

% user defined parallel operations must occur inside an spmd block
spmd
    fprintf(['hello! I am ' num2str(labindex) ' of ' num2str(numlabs) '\n']);
end


opt = odeset('reltol',1e-12,'abstol',1e-12);

r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
r2=[6729.43097302094	-318.159560024553	-1381.89229666774	2.26490482713380	-1.18931778372278	7.16940625613625]';


T1=[linspace(0,160690,10),[160691:160790]];
T2=[linspace(0,21550,10),[21551:21650]];

P1=blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
P2 =blkdiag(0.01,0.01,0.01,1e-8,1e-8,1e-8);

%%
if strcmp(method,'MC')
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
end


%%
if strcmp(method,'GH5')
    
    [X01gh5,w1gh5]=GH_points(r1,P1,5);
    [X02gh5,w2gh5]=GH_points(r2,P2,5);
    N=length(w1gh5);
    X1gh5=zeros(length(T1),6,N);
    X2gh5=zeros(length(T2),6,N);
    parfor i=1:N
        i
        [t,X1]= ode45(@twoBody,T1,X01gh5(i,:)',opt);
        [t,X2]= ode45(@twoBody,T2,X02gh5(i,:)',opt);
        X1gh5(:,:,i)=X1;
        X2gh5(:,:,i)=X2;
    end
    
    save('CollpaperSimsGH5','X1gh5','X2gh5','w1gh5','w2gh5','r1','r2','P1','P2','T1','T2','-v7.3')
    
end


%%
if strcmp(method,'GH6')
    
    [X01gh6,w1gh6]=GH_points(r1,P1,6);
    [X02gh6,w2gh6]=GH_points(r2,P2,6);
    N=length(w1gh6);
    X1gh6=zeros(length(T1),6,N);
    X2gh6=zeros(length(T2),6,N);
    parfor i=1:N
        i
        [t,X1]= ode45(@twoBody,T1,X01gh6(i,:)',opt);
        [t,X2]= ode45(@twoBody,T2,X02gh6(i,:)',opt);
        X1gh6(:,:,i)=X1;
        X2gh6(:,:,i)=X2;
    end
    
    save('CollpaperSimsGH6','X1gh6','X2gh6','w1gh6','w2gh6','r1','r2','P1','P2','T1','T2','-v7.3')
    
end

%%
if strcmp(method,'GH7')
    
    [X01gh7,w1gh7]=GH_points(r1,P1,7);
    [X02gh7,w2gh7]=GH_points(r2,P2,7);
    N=length(w1gh7);
    X1gh7=zeros(length(T1),6,N);
    X2gh7=zeros(length(T2),6,N);
    parfor i=1:N
        i
        [t,X1]= ode45(@twoBody,T1,X01gh7(i,:)',opt);
        [t,X2]= ode45(@twoBody,T2,X02gh7(i,:)',opt);
        X1gh7(:,:,i)=X1;
        X2gh7(:,:,i)=X2;
    end
    
    save('CollpaperSimsGH7','X1gh7','X2gh7','w1gh7','w2gh7','r1','r2','P1','P2','T1','T2','-v7.3')
    
end
%%
if strcmp(method,'CUT6')
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
    
     save('CollpaperSimsCUT6','X1cut6','X2cut6','w1cut6','w2cut6','r1','r2','P1','P2','T1','T2','-v7.3')
     
end
%%
if strcmp(method,'CUT8')
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
    save('CollpaperSimsCUT8','X1cut8','w1cut8','X2cut8','w2cut8','r1','r2','P1','P2','T1','T2','-v7.3')
    
end
%%
% all done, close pool of workers
% matlabpool close
% exit;