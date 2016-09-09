function MeasPairs=GetMeasPairs_satprob(MeasPairs,XsigSat,Radmodel,T0,Tk,Tvec,method,algo)
% tasking from T0 to Tk , including both T0 and Tk
% Assuminig T0-1 is already updated and current time step

% NO OPTIMIZATION ----------  EXHAUSTIVE SEARCH



global kappa
kappa=1;
switch lower(method)
    case 'ut'
        qd_pts=@(m,P)UT_sigmapoints(m,P,2);
    case 'cut4'
        qd_pts=@conjugate_dir_gausspts;
    case 'cut6'
        qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
    case 'cut8'
        qd_pts=@conjugate_dir_gausspts_till_8moment;
    case 'gh'
        qd_pts=@(m,P)GH_pts(m,P,para);
    otherwise
        keyboard
        error('smthg is wrong: DONT ask me what')
end

Nsat=Radmodel.Nsat;
Nrad=Radmodel.Nrad;
Tsim=Tvec(T0:Tk);


% update to T0 from T0-1
[XsigSat,~,~]=Propagation_Mu_Cov_satprob(1:Nsat,XsigSat,Radmodel,T0-1,T0,Tvec,method,{'NoProcessNoise'});

%% Greedy in time algorithm
% ATTENTION:  only works for 1 constraint problems on the sensor i.e.
% sensor can see only 1 target at a time

% [satellites,sensors/radars] :-- MeasPairs{Tsim(1)}

% Tasking starts from T0
% first greedy in time, then greedy in sensors
if max(strcmpi(algo,'GreedyTimeGreedySensors'))==1
    % greedy in time
    for tt=1:1:length(Tsim)
        
        MeasPairs{T0+tt-1}=[];
        
        % greedy in sensor
        for Nr=1:1:Nrad
            MI=zeros(1,Nsat);
            for Ns=1:1:Nsat
                % check visibility and jump if not a visible pair
                if check_rad_sat_pair_visibility(Ns,Nr,XsigSat,Radmodel,T0+tt-1,Tvec,method)<0
                    MI(Ns)=0;
                    continue
                end
                
                if isempty(MeasPairs{T0+tt-1})
                    k=[];
                else
                    k=find(MeasPairs{T0+tt-1}(:,1)==Ns); %[satellites,sensors/radars]
                end
                if isempty(k)
                    Nrads=Nr;
                else
                    Nrads=[MeasPairs{T0+tt-1}(k,2)',Nr]; %[satellites,sensors/radars]
                end
                tempMeasPairs=MeasPairs;
                
                tempMeasPairs{T0+tt-1}=[Ns*ones(length(Nrads),1),Nrads(:)];

                [~,~,Pk_t,~,Pku_t]=Meas_Update_satprob(Ns,XsigSat,tempMeasPairs,Radmodel,T0+tt-1,Tvec,method,[],{'pseudoupdate'});

                %                 det(Pk_t{Ns})/det(Pku_t{Ns})
                R=detratio(Pk_t{Ns},Pku_t{Ns});
                if R<1.01  %Pf/Pu
                    MI(Ns)=0;
                else
                    MI(Ns)=0.5*log(R);
                end
                %                 MI(Ns)=GetMI(Ns,Nrads,tt,Mu_t,P_t,Radmodel,qd_pts);
            end
            
            [~,ind]=max(MI);
            if isnan(MI(ind))==0 && MI(ind)>0
                
                MeasPairs{T0+tt-1}=vertcat(MeasPairs{T0+tt-1},[ind,Nr]);
            end
%             MeasPairs{T0+tt-1}
%             keyboard
        end
        
        % update the ones tasked and then propogate with onne time step using the measpairs sought
        [XsigSat,~,~,~,~]=Meas_Update_satprob(1:Nsat,XsigSat,MeasPairs,Radmodel,T0+tt-1,Tvec,method,[],{'pseudoupdate'});
        [XsigSat,~,~]=Propagation_Mu_Cov_satprob(1:Nsat,XsigSat,Radmodel,T0+tt-1,T0+tt,Tvec,method,{'NoProcessNoise'});
        
        
    end
    
    
end

end

% you can give multiple Nrads
function MI=GetMI(Nsat,Nrads,T_t,Mu_t,P_t,Radmodel,quadfunc)
Nrads=Nrads(:)';

mux=Mu_t{Nsat}(T_t,:)';
Px=reshape(P_t{Nsat}(T_t,:),6,6);
[X,w]=quadfunc(mux,Px);
N=size(X,1);

Y=[];
for nr=Nrads
    ZZ=zeros(N,Radmodel.hn);
    G=zeros(N,1);
    H=zeros(N,1);
    
    for msi=1:1:N
        ZZ(msi,:)=Radmodel.h(X(msi,:)',nr);
        [gg,hh]=Radmodel.G(X(msi,:)',nr);
        G(msi)=gg;
        H(msi)=hh;
    end
    
    if sum(isnan(H))<length(H)/2
        Y=horzcat(Y,ZZ);
        %                 RG=0;
        %                 for ii=1:1:N
        %                     RG=RG+w(ii)*G(ii)*Radmodel.R(Srad(nr));
        %                 end
        RG=Radmodel.R(nr);
        RR= blkdiag(RR,RG);
    end
end

if isempty(Y)
    MI=0;
else
    [mz,Pz]=MeanCov(Y,w);
    Pz=Pz+RR;
    Pcc=CrossCov(X,mux,Y,mz,w);
    K=Pcc/Pz;
    Pk1=Pxk-K*Pz*K';
    if det(Px)/det(Pk1)<1.01  %Pf/Pu
        MI=0;
    else
        MI=0.5*log(det(Px)/det(Pk1));
    end
end

end






