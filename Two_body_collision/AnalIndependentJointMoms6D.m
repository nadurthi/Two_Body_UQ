function [Y,MJ]=AnalIndependentJointMoms6D(Mx,My,N)

Dt=[ones(2,6)];
W=0.5*ones(2,1);
[y16D,M1tDth]=Cal_moments_samples(Dt ,W,1,'raw');
 [y26D,M2tDth]=Cal_moments_samples(Dt ,W,2,'raw');
 [y36D,M3tDth]=Cal_moments_samples(Dt ,W,3,'raw');
 [y46D,M4tDth]=Cal_moments_samples(Dt ,W,4,'raw');
 [y56D,M5tDth]=Cal_moments_samples(Dt ,W,5,'raw');
 [y66D,M6tDth]=Cal_moments_samples(Dt ,W,6,'raw');




if N>=1
MJ.M1=zeros(size(y16D,1),2);
Y.y1=y16D;
for i=1:1:size(y16D)
    ee=y16D(i,:);
    if sum(ee(1:3))==0
        MJ.M1(i,1)=1;
    else
        ind=MomentVecorder(ee(1:3));
        switch ind(1)
            case 1
                MJ.M1(i,1)=Mx.M1(ind(2));
            case 2
                MJ.M1(i,1)=Mx.M2(ind(2));
            case 3
                MJ.M1(i,1)=Mx.M3(ind(2));
            case 4
                MJ.M1(i,1)=Mx.M4(ind(2));
            case 5
                MJ.M1(i,1)=Mx.M5(ind(2));
            case 6
                MJ.M1(i,1)=Mx.M6(ind(2));          
        end
    
    end
     if sum(ee(4:6))==0
        MJ.M1(i,2)=1;
     else
        ind=MomentVecorder(ee(4:6));
        switch ind(1)
            case 1
                MJ.M1(i,2)=My.M1(ind(2));
            case 2
                MJ.M1(i,2)=My.M2(ind(2));
            case 3
                MJ.M1(i,2)=My.M3(ind(2));
            case 4
                MJ.M1(i,2)=My.M4(ind(2));
            case 5
                MJ.M1(i,2)=My.M5(ind(2));
            case 6
                MJ.M1(i,2)=My.M6(ind(2));          
        end    
     end 
end
MJ.M1=prod(MJ.M1,2);
end
%%
if N>=2
    MJ.M2=zeros(size(y26D,1),2);
Y.y2=y26D;
for i=1:1:size(y26D)
    ee=y26D(i,:);
    if sum(ee(1:3))==0
        MJ.M2(i,1)=1;
    else
        ind=MomentVecorder(ee(1:3));
        switch ind(1)
            case 1
                MJ.M2(i,1)=Mx.M1(ind(2));
            case 2
                MJ.M2(i,1)=Mx.M2(ind(2));
            case 3
                MJ.M2(i,1)=Mx.M3(ind(2));
            case 4
                MJ.M2(i,1)=Mx.M4(ind(2));
            case 5
                MJ.M2(i,1)=Mx.M5(ind(2));
            case 6
                MJ.M2(i,1)=Mx.M6(ind(2));          
        end
    
    end
     if sum(ee(4:6))==0
        MJ.M2(i,2)=1;
     else
        ind=MomentVecorder(ee(4:6));
        switch ind(1)
            case 1
                MJ.M2(i,2)=My.M1(ind(2));
            case 2
                MJ.M2(i,2)=My.M2(ind(2));
            case 3
                MJ.M2(i,2)=My.M3(ind(2));
            case 4
                MJ.M2(i,2)=My.M4(ind(2));
            case 5
                MJ.M2(i,2)=My.M5(ind(2));
            case 6
                MJ.M2(i,2)=My.M6(ind(2));          
        end    
     end 
end
MJ.M2=prod(MJ.M2,2);
end
%%
if N>=3
    MJ.M3=zeros(size(y36D,1),2);
Y.y3=y36D;
for i=1:1:size(y36D)
    ee=y36D(i,:);
    if sum(ee(1:3))==0
        MJ.M3(i,1)=1;
    else
        ind=MomentVecorder(ee(1:3));
        switch ind(1)
            case 1
                MJ.M3(i,1)=Mx.M1(ind(2));
            case 2
                MJ.M3(i,1)=Mx.M2(ind(2));
            case 3
                MJ.M3(i,1)=Mx.M3(ind(2));
            case 4
                MJ.M3(i,1)=Mx.M4(ind(2));
            case 5
                MJ.M3(i,1)=Mx.M5(ind(2));
            case 6
                MJ.M3(i,1)=Mx.M6(ind(2));          
        end
    
    end
     if sum(ee(4:6))==0
        MJ.M3(i,2)=1;
     else
        ind=MomentVecorder(ee(4:6));
        switch ind(1)
            case 1
                MJ.M3(i,2)=My.M1(ind(2));
            case 2
                MJ.M3(i,2)=My.M2(ind(2));
            case 3
                MJ.M3(i,2)=My.M3(ind(2));
            case 4
                MJ.M3(i,2)=My.M4(ind(2));
            case 5
                MJ.M3(i,2)=My.M5(ind(2));
            case 6
                MJ.M3(i,2)=My.M6(ind(2));          
        end    
     end 
end
MJ.M3=prod(MJ.M3,2);
end
%%
if N>=4
    MJ.M4=zeros(size(y46D,1),2);
Y.y4=y46D;
for i=1:1:size(y46D)
    ee=y46D(i,:);
    if sum(ee(1:3))==0
        MJ.M4(i,1)=1;
    else
        ind=MomentVecorder(ee(1:3));
        switch ind(1)
            case 1
                MJ.M4(i,1)=Mx.M1(ind(2));
            case 2
                MJ.M4(i,1)=Mx.M2(ind(2));
            case 3
                MJ.M4(i,1)=Mx.M3(ind(2));
            case 4
                MJ.M4(i,1)=Mx.M4(ind(2));
            case 5
                MJ.M4(i,1)=Mx.M5(ind(2));
            case 6
                MJ.M4(i,1)=Mx.M6(ind(2));          
        end
    
    end
     if sum(ee(4:6))==0
        MJ.M4(i,2)=1;
     else
        ind=MomentVecorder(ee(4:6));
        switch ind(1)
            case 1
                MJ.M4(i,2)=My.M1(ind(2));
            case 2
                MJ.M4(i,2)=My.M2(ind(2));
            case 3
                MJ.M4(i,2)=My.M3(ind(2));
            case 4
                MJ.M4(i,2)=My.M4(ind(2));
            case 5
                MJ.M4(i,2)=My.M5(ind(2));
            case 6
                MJ.M4(i,2)=My.M6(ind(2));          
        end    
     end 
end
MJ.M4=prod(MJ.M4,2);
end
%%
if N>=5
    MJ.M5=zeros(size(y56D,1),2);
Y.y5=y56D;
for i=1:1:size(y56D)
    ee=y56D(i,:);
    if sum(ee(1:3))==0
        MJ.M5(i,1)=1;
    else
        ind=MomentVecorder(ee(1:3));
        switch ind(1)
            case 1
                MJ.M5(i,1)=Mx.M1(ind(2));
            case 2
                MJ.M5(i,1)=Mx.M2(ind(2));
            case 3
                MJ.M5(i,1)=Mx.M3(ind(2));
            case 4
                MJ.M5(i,1)=Mx.M4(ind(2));
            case 5
                MJ.M5(i,1)=Mx.M5(ind(2));
            case 6
                MJ.M5(i,1)=Mx.M6(ind(2));          
        end
    
    end
     if sum(ee(4:6))==0
        MJ.M5(i,2)=1;
     else
        ind=MomentVecorder(ee(4:6));
        switch ind(1)
            case 1
                MJ.M5(i,2)=My.M1(ind(2));
            case 2
                MJ.M5(i,2)=My.M2(ind(2));
            case 3
                MJ.M5(i,2)=My.M3(ind(2));
            case 4
                MJ.M5(i,2)=My.M4(ind(2));
            case 5
                MJ.M5(i,2)=My.M5(ind(2));
            case 6
                MJ.M5(i,2)=My.M6(ind(2));          
        end    
     end 
end
MJ.M5=prod(MJ.M5,2);
end
%%
if N>=6
    MJ.M6=zeros(size(y66D,1),2);
    Y.y6=y66D;
for i=1:1:size(y66D)
    ee=y66D(i,:);
    if sum(ee(1:3))==0
        MJ.M6(i,1)=1;
    else
        ind=MomentVecorder(ee(1:3));
        switch ind(1)
            case 1
                MJ.M6(i,1)=Mx.M1(ind(2));
            case 2
                MJ.M6(i,1)=Mx.M2(ind(2));
            case 3
                MJ.M6(i,1)=Mx.M3(ind(2));
            case 4
                MJ.M6(i,1)=Mx.M4(ind(2));
            case 5
                MJ.M6(i,1)=Mx.M5(ind(2));
            case 6
                MJ.M6(i,1)=Mx.M6(ind(2));          
        end
    
    end
     if sum(ee(4:6))==0
        MJ.M6(i,2)=1;
     else
        ind=MomentVecorder(ee(4:6));
        switch ind(1)
            case 1
                MJ.M6(i,2)=My.M1(ind(2));
            case 2
                MJ.M6(i,2)=My.M2(ind(2));
            case 3
                MJ.M6(i,2)=My.M3(ind(2));
            case 4
                MJ.M6(i,2)=My.M4(ind(2));
            case 5
                MJ.M6(i,2)=My.M5(ind(2));
            case 6
                MJ.M6(i,2)=My.M6(ind(2));          
        end    
     end 
end
MJ.M6=prod(MJ.M6,2);
end