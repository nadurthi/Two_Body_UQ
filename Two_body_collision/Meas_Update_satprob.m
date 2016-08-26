function [XsigSat,Muk_t,Pk_t,Muku_t,Pku_t]=Meas_Update_satprob(Nsats,XsigSat,MeasPairs,Radmodel,Tk,Tvec,method,ytruth,further)
% Nsats is the index of the satellites to be updated
% simply do the measurement update for all the satellites at the time step
% Tk

% make the measurement if ture, or just use the pseudo measurements
% Tk is absolute time in the Tvec

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
        error('smthg is wrong: DONT ask me what')
end


Nsat=Radmodel.Nsat;
Tsim=Tvec(Tk);

Muk_t=cell(Nsat,1);
Pk_t=cell(Nsat,1);
Muku_t=cell(Nsat,1);
Pku_t=cell(Nsat,1);

if max(strcmpi(further,'TaskAll'))==1
    % overite teh tasking to task all sensors to all the satellites
    MeasPairs{Tk}=[];
    for Ns=Nsats
        MeasPairs{Tk}=vertcat(MeasPairs{Tk},[Ns*ones(Radmodel.Nrad,1) ,[1:Radmodel.Nrad]']);
    end
end


for Ns=Nsats
    
    mk=XsigSat{Ns,1}(Tk,:)';
    Pk=reshape(XsigSat{Ns,2}(Tk,:),6,6);
    
    Muk_t{Ns}=mk;
    Pk_t{Ns}=Pk;
    
    [x,w]=qd_pts(mk,Pk);
    N=size(x,1);
    
    if sum(size(MeasPairs{Tk}))==0
        ind=[];
    else
        ind=find(MeasPairs{Tk}(:,1)==Ns);  %[satellites,sensors/radars]
    end
    
    
    
    if isempty(ind)==0
        Srad=MeasPairs{Tk}(ind,2);
        %         Y=zeros(size(x,1),Radmodel.hn*length(Srad));
        Y=[];
        ym=[];
        RR=[];
        GG=[];
        for nr=1:1:length(Srad)
            ZZ=zeros(N,Radmodel.hn);
            G=zeros(N,1);
            H=zeros(N,1);
            for msi=1:1:N
                if isreal(x(msi,:))==0
                    keyboard
                end
%                 try
                ZZ(msi,:)=Radmodel.h(x(msi,:)',Srad(nr));
%                 catch
%                     keyboard
%                 end
                   
                [gg,hh]=Radmodel.G(x(msi,:)',Srad(nr));
                G(msi)=gg;
                H(msi)=hh;
                %                 if isnan(ZZ(msi,1))==1
                %                     flag1=1;
                %                     break;
                %                 end
            end
            %             if flag1==0
            %             Y(msi,(nr-1)*Radmodel.hn+1:(nr*Radmodel.hn))=ZZ(:)';

            if sum(isnan(H))<length(H)/2
                Y=horzcat(Y,ZZ);
                RG=0;
                %                 for ii=1:1:N
                %                     RG=RG+w(ii)*G(ii)*Radmodel.R(Srad(nr));
                %                 end
                RG=Radmodel.R(Srad(nr));
                RR= blkdiag(RR,RG);
                if max(strcmpi(further,'pseudoupdate'))==1
                    ym=-1234.1234;
                else
                    yjk=Radmodel.h(ytruth{Ns,1}(Tk,:),Srad(nr))+sqrtm(Radmodel.R(Srad(nr)))*randn(Radmodel.hn,1);
                    ym=vertcat(ym,yjk);
                end
            end
        end
        
        
        
        
        if isempty(ym)==0 && sum(isnan(ym))==0
            [mz,Pz]=MeanCov(Y,w);
            Pz=Pz+RR;
            Pcc=CrossCov(x,mk,Y,mz,w);
            if max(strcmpi(further,'pseudoupdate'))==1
                K=Pcc/Pz;
                Pku=Pk-K*Pz*K';
                mku=mk;
                
            else
                disp(strcat('Measurement Update for sat ',num2str(Ns)))
                [mku,Pku]=KalmanUpdate(mk,Pk,mz,Pz,Pcc,ym);
            end
            
            if max(strcmpi(further,'NoLoad'))==1
                
            else
                XsigSat{Ns,1}(Tk,:)=mku';
                XsigSat{Ns,2}(Tk,:)=reshape(Pku,1,36);
            end
            Muku_t{Ns}=mku;
            Pku_t{Ns}=Pku;
            
        else
            Muku_t{Ns}=mk;
            Pku_t{Ns}=Pk;
        end
        
    else
        Muku_t{Ns}=mk;
        Pku_t{Ns}=Pk;
        
    end
    
    
end
