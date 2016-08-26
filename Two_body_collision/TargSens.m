%% A class that containts all targets and all sensors
%         SensMsgs={'OutOfFOV','PseudoUpdate','SimpleTask'};
%         SensStates={'HasToTakeMeas','AlreadyTakenMeas'}
%         SensTypes={'Move','Stationary'};
%         TargTypes={'Move','Stationary','Virtual'};
%         TargStates={'CurrAt_tk','aPrioriAt_tk','aPostAt_tk'};
%         UQType={'UT','CUT4','CUT6','CUT8'};
%         InfoType={'FIM','MI'};
%         AvailSensors={'Range+Bearing','Range','Bearing'};
%         AvailSensors_dim={2,1,1};
%         AvailTargDyn={'UM','CT','NoDyn'};
%         AvailTargDyn_dim={4,5,2};
%         FOVpenaltytypes={'Simple','QuadraticOUT','NoPenalty'};
%         SensConstraintTypes={'1Sensor->1Target','1Sensor->AllTarget','1Sensor->1Target/AllTarget'};
%
%        If a sensor has more than 1 mode , add it as another sensor. Then
%        in the tasking phase put a constraint that only one of these
%        sensors can operate at a given time. The sensor variables are
%        arranged as
%              T        T       T      C        C    ||
%         --   s1      s2      s3      s4      s5    ||
%         t1   s11     s21     s31     s41     s51   ||  1    5   9    13 .
%         t2   s12     s22     s32     s42     s52   ||  2    6   10   14 .
%         t3   s13     s23     s33     s43     s53   ||  3    7   11   15 .
%         t4   s14     s24     s34     s44     s54   ||  4    8   12   16 .

% T means a traking sensors that can see only one object at once
% C means a coverage sensor ... i.e. can see all within the FOV
% Constraints: s11+s12+s13+s14<=1  for T sensors
% S41,s42,s43,s44   for eg: s43=0 if sensor 4 cannot see the 3rd target
% The variables are stacked as [s11,s12,s13,s14  ,s21,s22,s23,s24
% ,s31,s32,s33,s34  ,....] Nt x Ns

classdef TargSens
    properties
        xlim
        ylim
        Ns
        Nt
        Ntvec_deleted  %operate only on this vector of valid targets
        x0
        P0
        xk  %current state estimates of the targets    xk={{k}{targid}}
        Pk  % Pk={{k}{targid}}
        xktruth
        hn
        fn
        sk  % current position of the sensor
        alphak % half anfle of FOV
        rmaxk  % max range of FOV
        dirk  % pointing angle of sensor
        TargType   % (obj.Nt,[1,2])=(Move/Station,Dynamic Model) indexes
        TargState  %{time}(time,:)current state : is is propogated, is it meas updated, or nothing done,
        % this isvec of structs
        SensState   %current state : has to take measurement, already taken meas...
        SensType  % {Move/Station,Sensor Model}
        SensConstraints  % description of what kind of sensor it is
        SensModeConstraintsVec  %if sensor is multimodal  [id of sens1_mode1,id of sens1_mode2; ... ] so mode1 and mode2 cannot work together
        Nmodeconstraints % just the size of SensModeConstraintsVec
        MU  %sensor-target decision  matrix  targs along rows, sesnsor along columns
        MUvec % the vec version of MU decision matrix
        MUindex  % matrix that contains the indices of the MUvec
        R
        Q
        G  %penalty for FOV out sens
        t0
        tk  % the current tiems step
        tf
        nt
        Nt_t0  %the starting time of the target. esp for targets added in the middle
        Tvec
        dt
        SensTargtask  %[sensors,targets]
        yk % current measurement
        quadpts
        A  %As<=b constraints
        b
        C %Cs==d
        d
    end
    properties(Constant)
        SensMsgs={'OutOfFOV','PseudoUpdate','SimpleTask'};
        
        SensStates={'HasToTakeMeas','AlreadyTakenMeas'}
        SensTypes={'Move','Stationary'};
        TargTypes={'Move','Stationary','Virtual'};
        TargStates={'CurrAt_tk','aPrioriAt_tk','aPostAt_tk'};
        
        UQType={'UT','CUT4','CUT6','CUT8'};
        InfoType={'FIM','MI'};
        AvailSensors={'Range+Bearing','Range','Bearing'};
        AvailSensors_dim={2,1,1};
        AvailTargDyn={'UM','CT','NoDyn'};
        AvailTargDyn_dim={4,5,2};
        
        FOVpenaltytypes={'Simple','QuadraticOUT','NoPenalty'};
        
        SensorTaskingMethods={'MIUB','FIM'}
        
        SensConstraintTypes={'1Sensor->1Target','1Sensor->AllTarget','1Sensor->AllTarget/1Target'};
        % 1 sensor can see only 1 target within FOV
        % 1 sensor can see all target within FOV
        % 1 sensor can choose see only 1 target or all target  within FOV
        
    end
    methods
        
        %============================================================================
        function obj=TargSens(t0,dt,tf)
            obj.Ns=0;
            obj.Nt=0;
            obj.t0=t0;
            obj.tf=tf;
            obj.dt=dt;
            obj.Tvec=t0:dt:tf;
            obj.nt=length(obj.Tvec);
            obj.tk=1;  %set the curr time to 0 or Tvec(1)
            obj.x0={};
            obj.P0={};
            obj.xk=cell(1,1); %xk{obj.Nt}=[xk0;xk1;xk2].....
            obj.Pk=cell(1,1);
            obj.xktruth=cell(1,1);
            obj.TargType=[0,0];  %(Move/Station,Dynamic Model)
            obj.SensType=[0,0];  %(Move/Station,Sensor Model)
            
            obj.TargState=cell(1,1);
            obj.SensState=cell(1,1);   %current state : has to take measurement, already taken meas...
            
            obj.R=cell(1,1); %one for each sensor
            obj.Q=cell(1,1);  % one for each target
            obj.yk=cell(1,1);  % each time the measurements are loaded into this
            obj.SensTargtask=cell(1,1); %each cell  [sensors,targets] is a vector
            obj.hn=0;
            obj.fn=0;
            obj.Ntvec_deleted=[];
            obj.Nt_t0=0;
            obj.sk=cell(1,1);
            obj.alphak=cell(1,1);
            obj.rmaxk=cell(1,1);
            obj.dirk=cell(1,1);
            
            obj.SensConstraints=cell(1,1);
            obj.C=[];
            obj.A=[];
            obj.b=[];
            obj.d=[];
            
            obj.SensModeConstraintsVec=cell(1,1);
            obj.Nmodeconstraints=0;
            obj.MUindex=0;
        end
        %============================================================================
        function index=Index_MUvec(obj,targid,sensid)
            % get the index number
            if size(obj.MUindex,1)~=obj.Nt || size(obj.MUindex,2)~=obj.Ns % checking if this index is already created or not
                obj.MUindex=zeros(obj.Nt,obj.Ns);
                k=1;
                for j=1:1:obj.Ns
                    for i=1:1:obj.Nt
                        obj.MUindex(i,j)=k;
                        k=k+1;
                    end
                end
            end
            if isempty(targid) %give all the indexes of the sensor sensid
                index=obj.MUindex(:,sensid);
                return
            end
            if isempty(sensid) %give all the indexes of the sensor sensid
                index=obj.MUindex(targid,:);
                return
            end
            index=obj.MUindex(targid,sensid);
            
            
        end
        %============================================================================
        function obj=Set_Sigmapts(obj,type)
            
            switch(find(strcmpi(obj.UQType,type)==1))
                case 1   %UT
                    obj.quadpts=@(mu,P)UT_sigmapoints(mu(:),P,2);
                case 2   %CUT4
            end
        end
        %============================================================================
        function obj=Add_Target(obj,TargType,DynType,xk,Pk,truth,Qk)  % targtype=Move/Stationary , truth is for all times
            %this way the target can be initiated at any time
            if max(strcmpi(TargType,obj.TargTypes))==0
                error('Wrong Targtype');
            end
            if max(strcmpi(DynType,obj.AvailTargDyn))==0
                error('Wrong DynType');
            end
            
            obj.Nt=obj.Nt+1;
            
            
            obj.TargType(obj.Nt,1)=find(strcmpi( obj.TargTypes,TargType)==1);
            obj.TargType(obj.Nt,2)=find(strcmpi( obj.AvailTargDyn,DynType)==1);
            obj.fn(obj.Nt)=obj.AvailTargDyn_dim{find(strcmpi( obj.AvailTargDyn,DynType)==1)};
            
            obj.TargState{obj.Nt}(obj.tk)=1;  % 1 is 'CurrAt_tk'
            
            if obj.tk==1
                obj.x0{obj.Nt}=xk(:)';
                obj.P0{obj.Nt}=reshape(Pk(:),1,obj.fn(obj.Nt)^2);
            end
            obj.xk{obj.Nt}(obj.tk,:)=xk(:)';
            obj.Pk{obj.Nt}(obj.tk,:)=reshape(Pk(:),1,obj.fn(obj.Nt)^2);
            obj.Q{obj.Nt}=Qk;
            
            obj.xktruth{obj.Nt}=truth; % for all times has to be given
            if size(truth,1)<obj.nt
                error('truth has to be given for all times')
            end
        end
        %============================================================================
        function obj=Delete_Target(obj,targid)  % list of deleted targets
            obj.Ntvec_deleted=horzcaT(obj.Ntvec_deleted,targid);
            
        end
        %============================================================================
        function obj=Add_Sensor(obj,SensorType,SensorModel,SensorConstraints,R,sk,alphak,rmaxk,dirn,FOVpenaltytype)
            
            
            if max(strcmpi(SensorType,obj.SensTypes))==0
                error('Wrong SensTypes');
            end
            if max(strcmpi(SensorModel,obj.AvailSensors))==0
                error('Wrong Sensor model');
            end
            
            obj.Ns=obj.Ns+1;
            
            obj.SensType(obj.Ns,1)=find(strcmpi( obj.SensTypes,SensorType)==1);%(Move/Station,Sensor Model)
            obj.SensType(obj.Ns,2)=find(strcmpi( obj.AvailSensors,SensorModel)==1);
            obj.hn(obj.Ns)=obj.AvailSensors_dim{find(strcmpi( obj.AvailSensors,SensorModel)==1)};
            
%                keyboard
            obj.sk{obj.Ns}(obj.tk,:)=sk;
            obj.alphak{obj.Ns}(obj.tk)=alphak;
            obj.rmaxk{obj.Ns}(obj.tk)=rmaxk;
            obj.dirk{obj.Ns}(obj.tk)=dirn;
            
            obj.SensState{obj.Ns}(obj.tk)=1; % has to make measurement
            
            obj.R{obj.Ns}=R;
            obj.G(obj.Ns)=find(strcmpi( obj.FOVpenaltytypes,FOVpenaltytype)==1);
            
            obj.SensConstraints{obj.Ns}=find(strcmpi( obj.SensConstraintTypes,SensorConstraints)==1);
            
        end
        %============================================================================
        function obj=Add_MultiModalSensor(obj,nmodes,SensorType,SensorModel,SensorConstraints,R,sk,alphak,rmaxk,dirn,FOVpenaltytype)
            % add all the modes of this sensor
            % give the properties in a cell in an order
            V=[];
            for nm=1:1:nmodes
             
                obj=Add_Sensor(obj,SensorType{nm},SensorModel{nm},SensorConstraints{nm},R{nm},sk,alphak{nm},rmaxk{nm},dirn{nm},FOVpenaltytype);
                V=horzcat(V,obj.Ns);
            end
            
            obj.Nmodeconstraints=obj.Nmodeconstraints+1;
            obj.SensModeConstraintsVec{obj.Nmodeconstraints}=V;
            %             coeff=zeros(1,obj.Ns);
            %             coeff(end-nmodes:1:end)=1;
            %             %Add the contraint that only mode can work
            %             obj=Add_InterSensConstraints(obj,0,1,coeff,'<=');
        end
        %============================================================================
        function [y,H,G]=SensorMeasurement_pts(obj,Tk,X,sensid)
            % for all the points in X at time Tk
            alpha=obj.alphak{sensid}(Tk);
            Rmax=obj.rmaxk{sensid}(Tk);
            dirn=obj.dirk{sensid}(Tk);
            xsenspos=obj.sk{sensid}(Tk,:);
            
            y=zeros(size(X,1),obj.hn(sensid));
            G=zeros(size(X,1),1);
            H=ones(size(X,1),1); % ones means inside
            
            for i=1:1:size(X,1)
                
                r=norm([X(i,1),X(i,2)]-[xsenspos(1),xsenspos(2)]);
                th=atan2(X(i,2)-xsenspos(2),X(i,1)-xsenspos(1));
                %alpha is (-pi,pi)
                %dirn is (-pi,pi)
                %th is (-pi,pi)
                
                diff=dirn-th;
                if diff>pi
                    diff=diff-2*pi;
                end
                if diff<-pi
                    diff=diff+2*pi;
                end
                
                switch obj.SensType(sensid,2)
                    case 1   %range and bearing
                        y(i,:)=[r,th];
                    case 2   % only range
                        y(i,:)=r;
                    case 3   % only bearing
                        y(i,:)=th;
                    otherwise
                        error('Cannot find sensor model')
                end
                
                G(i)=0;
                if r>Rmax
                    switch obj.G(obj.Ns)
                        case 1   % simple constant penalty
                            G(i)=10;
                        case 2   % quadratic penalty
                            G(i)=5*(r-Rmax)^2;
                        case 3   % no penalty
                            G(i)=1;
                    end
                    H(i)=0;  % 0 means the reading is outside the FOV
                end
                if abs(diff)>alpha
                    switch obj.G(obj.Ns)
                        case 1   % simple constant penalty
                            G(i)=G(i)+10;
                        case 2   % quadratic penalty
                            G(i)=G(i)+5*(abs(diff)-alpha)^2;
                        case 3   % no penalty
                            G(i)=1;
                    end
                    H(i)=0;  % 0 means the reading is outside the FOV
                end
                if abs(diff)<alpha && r<Rmax
                    G(i)=1;
                end
            end
            
        end
        %============================================================================
        function obj=DynProp_Sensor(obj,Tk,sensid)
            % update to-> Tk+1 from Tk values without modifying the tk value
            sensid=sensid(:)';
            for ss=sensid
                switch obj.SensType(ss,1)  % first element is the move/stationary
                    case 1   % the special case where sensors move
                        error('Method of dyn sens not implemented yet');
                    case 2   %case 2 is stationary
                        obj.alphak{ss}(Tk+1)=obj.alphak{ss}(Tk);
                        obj.rmaxk{ss}(Tk+1)=obj.rmaxk{ss}(Tk);
                        obj.dirk{ss}(Tk+1)=obj.dirk{ss}(Tk);
                        obj.sk{ss}(Tk+1,:)=obj.sk{ss}(Tk,:);
                end
            end
        end
        %============================================================================
        function obj=UpdateTimek(obj)
            obj.tk=obj.tk+1;
        end
        %============================================================================
        function obj=DynProp_Target(obj,targs)
            % prop grom tk to tk+1
            %the time variable is not updated
            
            %first remove the targets that have been deleted i.e. do
            %nothing to them
            for targid=targs(:)'
                if isempty(find(obj.Ntvec_deleted==targid))==0
                    continue;
                end
                %                 if obj.TargType(targid,1)==2  % 2 is stationary so just copy the value into the next
                %                     Xfk=obj.xk{targid}(obj.tk,:);
                %                     Pk1=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
                % %                     Pk1=Pk1+obj.Q{targid};
                % %                     obj.Pk{targid}(obj.tk+1,:)=reshape(Pk1,1,obj.fn(targid)^2);
                % %                     continue;
                %                 end
                
                Xfk=obj.xk{targid}(obj.tk,:)';
                Pfk=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
                
                [Xk,wk]=obj.quadpts(Xfk,Pfk);
                Xk1=zeros(size(Xk));
                for i=1:1:length(wk)
                    switch(obj.TargType(targid,2)) % type of dynamics
                        case 1  % UM motion
                            Xk1(i,:)=KIRB_UM_eg_dyn_disc(Xk(i,:)',obj.dt);
                        case 2  % CT motion
                            Xk1(i,:)=KIRB_CT_eg_dyn_disc(Xk(i,:)',obj.dt);
                        case 3  % No dynamics case
                            Xk1(i,:)=Xk(i,:);
                    end
                end
                
                [xk1,Pk1]=MeanCov(Xk1,wk);
                Pk1=Pk1+obj.Q{targid};
                
                obj.xk{targid}(obj.tk+1,:)=xk1;
                obj.Pk{targid}(obj.tk+1,:)=reshape(Pk1,1,obj.fn(targid)^2);
            end
        end
        %============================================================================
        function obj=MeasUpdate_Target(obj,targs)
            %             do the measurement update for specified targets at curr time
            %             tk, with out modifying the time variable
            % RMEMEBER SENS-TASKING HAS TO BE DONE
            % SensTargtask{tk}(sensors,target)
            
            for targid=targs(:)'
                %Do nothing to the deleted targets
                if isempty(find(obj.Ntvec_deleted==targid))==0
                    continue;
                end
                
                sensind=find(obj.SensTargtask{obj.tk}(:,2)==targid);
                if isempty(sensind)
                    continue;
                end
                Xfk=obj.xk{targid}(obj.tk,:)';
                Pfk=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
                [Xk,wk]=obj.quadpts(Xfk,Pfk);
                
                % Stack all the measurmeents together for the Kalman Update
                RR=[];
                YM=[];
                Y=[];
                sensids=obj.SensTargtask{obj.tk}(sensind,1);
                for nsens=sensids(:)'
                    [yy,H,G]=SensorMeasurement_pts(obj,obj.tk,Xk,nsens);
                    if length(find(H==0))==length(H)  %if no. of outofFOV sigpts == no. of sigpts
                        continue; % if no sig pt is visible just dont use this sensor
                    end
                    [ym,~,~]=SensorMeasurement_pts(obj,obj.tk,obj.xktruth{targid}(obj.tk,:),nsens);
                    ym=ym+sqrtm(obj.R{nsens})*randn(obj.hn(obj.Ns),1);
                    YM=vertcat(YM,ym);
                    Y=horzcat(Y,yy);
                    RR=blkdiag(RR,obj.R{nsens});
                end
                
                if isempty(Y)==0
                    [mz,Pz]=MeanCov(Y,wk);
                    Pz=Pz+RR;
                    Pcc=CrossCov(Xk,Xfk,Y,mz,wk);
                    [xku,Pku]=KalmanUpdate(Xfk,Pfk,mz,Pz,Pcc,YM);
                else
                    xku=Xfk;
                    Pku=Pfk;
                end
                obj.xk{targid}(obj.tk,:)=xku;
                obj.Pk{targid}(obj.tk,:)=reshape(Pku,1,obj.fn(targid)^2);
            end
            
        end
        %============================================================================
        function [Pfk,Pku, Xfk,xku]=MeasUpdate_Pseudo(obj,targid,sensid)
            %  Do simple pseudo measurement update to get the updated mean
            %SINGLE TARGET ONLY
            % and covarianr for tasking purposes
            % DO NOT CHANGE TIME OR UPDATE OBJ states
            if isempty(find(obj.Ntvec_deleted==targid))==0
                Pfk=0;
                Pku=0;
                Xfk=0;
                xku=0;
                disp(['target ',num2str(targid),' has been deleted nothing is being done to it'])
                return
            end
            
            Xfk=obj.xk{targid}(obj.tk,:)';
            Pfk=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
            try
            [Xk,wk]=obj.quadpts(Xfk,Pfk);
            catch
                keyboard
            end
            % Stack all the measurmeents together for the Kalman Update
            RR=[];
            Y=[];
            
            for nsens=sensid(:)'
                [yy,H,GG]=SensorMeasurement_pts(obj,obj.tk,Xk,nsens);
                if length(find(H==0))==length(H)  %if no. of outofFOV sigpts == no. of sigpts
%                     disp(['no sig pt of target ',num2str(targid),' is visible just dont use sensor ',num2str(sensid)])
                    continue; % if no sig pt is visible just dont use this sensor
                end
                
                Y=horzcat(Y,yy);
                Ri=0;
                for i=1:1:length(wk)
                    Ri=Ri+wk(i)*GG(i)^2*obj.R{nsens};
                end
                RR=blkdiag(RR,Ri);
            end
            
            if isempty(Y)==0
                [mz,Pz]=MeanCov(Y,wk);
                Pz=Pz+RR;
                Pcc=CrossCov(Xk,Xfk,Y,mz,wk);
                [xku,Pku]=KalmanUpdate(Xfk,Pfk,mz,Pz,Pcc,mz); % use mz as measurement to cancle the innovation
                 disp(['Updated target ',num2str(targid),' with sensor ',num2str(sensid)])
                 
            else
                xku=Xfk;
                Pku=Pfk;
               
            end
        end
        %============================================================================
        function obj=OptiSens_TargPairs(obj,type)
            %             type= 'Simple'... all possible .. no optimization but respect
            %             the FOV constraints: I think i Will implement
            %             this as a greedy algorithms
            %             type = 'optimized' Solve an optimization problem
            %             obj.SensTargtask{obj.tk}=[0 0] when there are no pairs atall
            %             ADD ALL THE REQUIRED CONSTRAINTS here
            %
            obj.SensTargtask{obj.tk}=[];
            obj.C=[];
            obj.A=[];
            obj.b=[];
            obj.d=[];
            II=eye(obj.Nt*obj.Ns);
            % adding constraints to deselect all the deleted target
            % variables
            for targid=1:1:obj.Nt
                if isempty(find(obj.Ntvec_deleted==targid))==0
                    ind=Index_MUvec(obj,targid,[]);
                    for i=1:1:length(ind)
                        c=zeros(1,obj.Ns*obj.Nt);
                        c(ind(i))=1;
                        obj.C=vertcat(obj.C,c);
                        obj.d=vertcat(obj.d,0);
                    end
                end
                
            end
            
            % adding the T and C constraints
            for sensid=1:1:obj.Ns
                if obj.SensConstraints{sensid}==1  % tracking sensor
                    ind=Index_MUvec(obj,[],sensid);
                    c=zeros(1,obj.Ns*obj.Nt);
                    c(ind)=1;
                    obj.A=vertcat(obj.A,c);
                    obj.b=vertcat(obj.b,1);
                end
                
            end
            % adding the mode cosntraints in obj.SensModeConstraintsVec
            % so that only 1 mode is selected
            for i=1:1:obj.Nmodeconstraints
                V=obj.SensModeConstraintsVec{i}; % V has the sensids that are related
                c=zeros(1,obj.Ns*obj.Nt);
                for j=1:1:length(V)
                    ind=Index_MUvec(obj,[],j);
                    if obj.SensConstraints{sensid}==1 % tracking sensor
                        c(ind)=1;
                    end
                    if obj.SensConstraints{sensid}==2 % coverage sensor
                        c(ind)=1/obj.Nt;
                    end
                end
                obj.A=vertcat(obj.A,c);
                obj.b=vertcat(obj.b,1);
            end
            %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if strcmpi(type,'MIUB') %implement this bintprog
                disp('MIub optimization is selected')
                % first for all pairs find the MI
                
                F=zeros(obj.Nt,obj.Ns);
                for targid=1:1:obj.Nt
                    for sensid=1:1:obj.Ns
                        F(targid,sensid)=MutualInformation(obj,targid,sensid);
                        if F(targid,sensid)==0  % remove this sensor by constraining to be 0
                            ind=Index_MUvec(obj,targid,sensid);
                            obj.C=vertcat(obj.C,II(ind,:));
                            obj.d=vertcat(obj.d,0);
                        end
                    end
                end
                % Next go thru all the coverage sensors and add that within
                % FOV variables are all equal
                for sensid=1:1:obj.Ns
                    if obj.SensConstraints{sensid}==2  % coverage sensor
                        F0=find(F(:,sensid)~=0);
                        if isempty(F0)
                            continue;
                        end
                        indMU=zeros(1,length(F0));
                        for i=1:1:length(F0)
                            indMU(i)=Index_MUvec(obj,F0(i),sensid);
                        end
                        c=zeros(length(F0),1);
                        c(1)=1; c(2)=-1;
                        Y=[];
                        for i=1:1:length(F0)
                            Y=horzcat(Y,circshift(c,i-1));
                        end
                        Y(end,:)=[];   % remove last cosntrraitn as s1=s2=s3 needs only 2 constraints
                        for i=1:1:length(F0)-1
                            c=zeros(1,obj.Ns*obj.Nt);
                            c(indMU)=Y(i,:);
                            obj.C=vertcat(obj.C,c);
                            obj.d=vertcat(obj.d,0);
                        end
                        
                    end
                    
                end
                keyboard
                %now solve the bintprog problem with the constraints
                f=reshape(F,obj.Ns*obj.Nt,1);
                obj.MUvec=bintprog(f,obj.A,obj.b,obj.C,obj.d);
                obj.MU=reshape(obj.MUvec,obj.Nt,obj.Ns);
                
                % now load the solutoin into the SensTask variabe;
                for targid=1:1:obj.Nt
                    for sensid=1:1:obj.Ns
                        if obj.MU(targid,sensid)==1
                            obj.SensTargtask{obj.tk}=vertcat(obj.SensTargtask{obj.tk},[sensid,targid]);   %[sensors,targets]
                        end
                    end
                end
                return
            end
            %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
        end
        %============================================================================
        function MI=MutualInformation(obj,targid,sensid)
            [Pf,Pu, ~,~]=MeasUpdate_Pseudo(obj,targid,sensid);
           
            if det(Pf)/det(Pu)<1.001
                MI=0;
            else
                MI=0.5*log(det(Pf)/det(Pu));
            end
        end
        %============================================================================
        function plot_scenario(obj) %No Tasking is considered
            figure(1)
            clf
            hold on
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.Ns    % all sensors are blue
                if obj.SensConstraints{sensid}==2  % it is a coverage sensor, plot circles
                    plot_ellipse(obj.sk{sensid}(obj.tk,1),obj.sk{sensid}(obj.tk,2),obj.rmaxk{sensid}(obj.tk),obj.rmaxk{sensid}(obj.tk),0,'b','fill');
                end
                if obj.SensConstraints{sensid}==1  % it is a track sensor, plot cones
                    plot_2Dcone(obj.sk{sensid}(obj.tk,1),obj.sk{sensid}(obj.tk,2),obj.alphak{sensid}(obj.tk),obj.rmaxk{sensid}(obj.tk),obj.dirk{sensid}(obj.tk),'g','fill')
                end
            end
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.Nt
                if isempty(find(obj.Ntvec_deleted==targid))==0
                    continue;
                end
                %plot covariance ellipse
                P=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
                plot_1sigellip(obj.xk{targid}(obj.tk,1:2),P(1:2,1:2),'r');
                %plot triangle as mean
                plot_triangle(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),2,1,0,'r','fill')
                %plot the trajectories est+truth
                plot(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),'r','linewidth',2)
                plot(obj.xktruth{targid}(obj.tk,1),obj.xktruth{targid}(obj.tk,2),'k--','linewidth',2)
            end
            
            
            hold off
        end
        %============================================================================
        function plot_Dynamic_scenario(obj) %Tasking is considered
            figure(1)
            clf
            hold on
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.Ns    % all sensors are blue
                if isempty(find(obj.SensTargtask{obj.tk}(:,1)==sensid))
                    continue;
                end
                if obj.SensConstraints{sensid}==2  % it is a coverage sensor, plot circles
                    plot_ellipse(obj.sk{sensid}(obj.tk,1),obj.sk{sensid}(obj.tk,2),obj.rmaxk{sensid}(obj.tk),obj.rmaxk{sensid}(obj.tk),0,'b','fill');
                end
                if obj.SensConstraints{sensid}==1  % it is a track sensor, plot cones
                    plot_2Dcone(obj.sk{sensid}(obj.tk,1),obj.sk{sensid}(obj.tk,2),obj.alphak{sensid}(obj.tk),obj.rmaxk{sensid}(obj.tk),obj.dirk{sensid}(obj.tk),'g','fill')
                end
            end
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.Nt
                if isempty(find(obj.Ntvec_deleted==targid))==0
                    continue;
                end
                %plot covariance ellipse
                P=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
                plot_1sigellip(obj.xk{targid}(obj.tk,1:2),P(1:2,1:2),'r');
                %plot triangle as mean
                plot_triangle(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),2,1,0,'r','fill')
                %plot the trajectories est+truth
                plot(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),'r','linewidth',2)
                plot(obj.xktruth{targid}(obj.tk,1),obj.xktruth{targid}(obj.tk,2),'k--','linewidth',2)
            end
            
            
            hold off
        end
        %============================================================================
        function plot_truth_scenario(obj) %No Tasking is considered+ ploy all truth
            figure(1)
            clf
            hold on
            % plotting the sensors and their FOV configurations
            for sensid=1:1:obj.Ns    % all sensors are blue
                if obj.SensConstraints{sensid}==2  % it is a coverage sensor, plot circles
                    
                    plot_ellipse(obj.sk{sensid}(obj.tk,1),obj.sk{sensid}(obj.tk,2),obj.rmaxk{sensid}(obj.tk),obj.rmaxk{sensid}(obj.tk),0,'b','fill');
                end
                if obj.SensConstraints{sensid}==1  % it is a track sensor, plot cones
                    plot_2Dcone(obj.sk{sensid}(obj.tk,1),obj.sk{sensid}(obj.tk,2),obj.alphak{sensid}(obj.tk),obj.rmaxk{sensid}(obj.tk),obj.dirk{sensid}(obj.tk),'g','fill')
                end
            end
            %plot the target with a triangle and an ellipse for confidence
            for targid=1:1:obj.Nt
                if isempty(find(obj.Ntvec_deleted==targid))==0
                    continue;
                end
                %plot covariance ellipse
                P=reshape(obj.Pk{targid}(obj.tk,:),obj.fn(targid),obj.fn(targid));
                plot_1sigellip(obj.xk{targid}(obj.tk,1:2),P(1:2,1:2),'r');
                %plot triangle as mean
                if obj.TargType(targid,1)==3  % if it is a vitual target
                    plot_triangle(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),1,1,0,'k','fill')
                    %plot the trajectories est+truth
                    plot(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),'r','linewidth',2)
                    plot(obj.xktruth{targid}(:,1),obj.xktruth{targid}(:,2),'k--','linewidth',2)
                else
                    plot_triangle(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),2,1,0,'r','fill')
                    %plot the trajectories est+truth
                    plot(obj.xk{targid}(obj.tk,1),obj.xk{targid}(obj.tk,2),'r','linewidth',2)
                    plot(obj.xktruth{targid}(:,1),obj.xktruth{targid}(:,2),'k--','linewidth',2)
                end
            end
            
%             axis([0,obj.xlim,0,obj.ylim])
            axis square
            axis equal
            hold off
            
        end
        %============================================================================
        function Detect_New_Targets(obj)
            % Make measurements only by coverage sensors only to find the
            % targets then add it to the list of targets
        end
        %============================================================================
    end
end