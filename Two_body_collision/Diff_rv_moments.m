function MDD=Diff_rv_moments(Mx,My)

d=length(Mx{1});
deg=size(Mx,1);

[X,w]=GH_points(zeros(d,1),eye(d),4);
[y1,M]=Cal_moments_samples(X,w,1,'raw');
[y2,M]=Cal_moments_samples(X,w,2,'raw');
[y3,M]=Cal_moments_samples(X,w,3,'raw');
[y4,M]=Cal_moments_samples(X,w,4,'raw');

MD=0;
lines=1;
if d==2
fid = fopen('diff_moms2D.m');
switch deg
    case 1
    totdeg=2;
    case 2
    totdeg=5;
    case 3
    totdeg=9;
    case 4
    totdeg=14;
end
end
if d==3
fid = fopen('diff_moms3D.m');
switch deg
    case 1
    totdeg=3;
    case 2
    totdeg=9;
    case 3
    totdeg=19;
    case 4
    totdeg=34;
end
end

dd=1;
while lines<=totdeg
    tline = fgets(fid);
    tline(strfind(tline,' '))=[];
    tline=strrep(tline, '+', ',+');
    tline=strrep(tline, '-', ',-');
    tline=strcat('+',tline);
        
    %splint the strings
        coms=strfind(tline,',');
        M=cell(length(coms),2);       
        kk=1;
        for j=1:1:length(coms)
           ss=tline(kk:coms(j)-1);
           if (strcmp(ss(2),'x')==1 || strcmp(ss(2),'y')==1)  
           M{j,1}=strcat(ss(1),'1');
           M{j,2}=ss(2:end);
           else
              M{j,1}=ss(1:2);
              M{j,2}=ss(3:end);
           end
           kk=coms(j)+1;
        end
        j=j+1;
        ss=tline(kk:end);
        if strcmp(ss(2),'x')==1 || strcmp(ss(2),'y')==1
           M{j,1}=strcat(ss(1),'1');
           M{j,2}=ss(2:end);
           else
              M{j,1}=ss(1:2);
              M{j,2}=ss(3:end);
        end
       %further split xs and ys as they are independent    
    for i=1:1:size(M,1)
    if isempty(strfind(M{i,2},'y'))
        if strcmp(M{i,2}(end-1),'^')==0
            M{i,2}=strcat(M{i,2},'^1');
        end
        continue;
    elseif isempty(strfind(M{i,2},'x'))
        if strcmp(M{i,2}(end-1),'^')==0
            M{i,2}=strcat(M{i,2},'^1');
        end
        M{i,3}=M{i,2};
        M{i,2}=[];
    else
        ind=strfind(M{i,2},'y');
        M{i,3}=M{i,2}(ind(1):end);
        M{i,2}(ind(1):end)=[];
        if strcmp(M{i,2}(end-1),'^')==0
            M{i,2}=strcat(M{i,2},'^1');
        end
        if strcmp(M{i,3}(end-1),'^')==0
            M{i,3}=strcat(M{i,3},'^1');
        end
    end
    end
    % now replacing them with their powers
    for i=1:1:size(M,1)
       D=[0,0,0];
       % x
        if isempty(M{i,2})==1
            M{i,2}=1; 
        else
            ind=strfind(M{i,2},'x1');
            if isempty(ind)
                D(1)=0;
            elseif strcmp(M{i,2}(ind+2),'^')
                 D(1)=str2num(M{i,2}(ind+3));
            else
                D(1)=1;
            end
            ind=strfind(M{i,2},'x2');
            if isempty(ind)
                D(2)=0;
            elseif strcmp(M{i,2}(ind+2),'^')
                 D(2)=str2num(M{i,2}(ind+3));
            else
                D(2)=1;
            end
            ind=strfind(M{i,2},'x3');
            if isempty(ind)
                D(3)=0;
            elseif strcmp(M{i,2}(ind+2),'^')
                 D(3)=str2num(M{i,2}(ind+3));
            else
                D(3)=1;
            end
            M{i,2}=D;
        end
        % y
        if isempty(M{i,3})==1
            M{i,3}=1; 
        else
            ind=strfind(M{i,3},'y1');
            if isempty(ind)
                D(1)=0;
            elseif strcmp(M{i,3}(ind+2),'^')
                 D(1)=str2num(M{i,3}(ind+3));
            else
                D(1)=1;
            end
            ind=strfind(M{i,3},'y2');
            if isempty(ind)
                D(2)=0;
            elseif strcmp(M{i,3}(ind+2),'^')
                 D(2)=str2num(M{i,3}(ind+3));
            else
                D(2)=1;
            end
            ind=strfind(M{i,3},'y3');
            if isempty(ind)
                D(3)=0;
            elseif strcmp(M{i,3}(ind+2),'^')
                 D(3)=str2num(M{i,3}(ind+3));
            else
                D(3)=1;
            end
            M{i,3}=D;
        end
      M{i,1}=str2num(M{i,1});
    end
    if d==2
    for i=1:1:size(M,1)
       if length(M{i,2})~=1
        M{i,2}(3)=[];
       end
       if length(M{i,3})~=1
        M{i,3}(3)=[];
       end
    end
    end
    % x moms
    for i=1:1:size(M,1)
       if length(M{i,2})==1
           continue;
       end
        switch sum(M{i,2})
           case 1
           M{i,2}=Mx{1}(find(sum(abs(repmat(M{i,2},size(y1,1),1)-y1),2)==0));    
            case 2
           M{i,2}=Mx{2}(find(sum(abs(repmat(M{i,2},size(y2,1),1)-y2),2)==0));  
            case 3
           M{i,2}=Mx{3}(find(sum(abs(repmat(M{i,2},size(y3,1),1)-y3),2)==0));
            case 4
           M{i,2}=Mx{4}(find(sum(abs(repmat(M{i,2},size(y4,1),1)-y4),2)==0));
       end
    end
        % y moms
    for i=1:1:size(M,1)
       if length(M{i,3})==1
           continue;
       end
        switch sum(M{i,3})
           case 1
           M{i,3}=My{1}(find(sum(abs(repmat(M{i,3},size(y1,1),1)-y1),2)==0));    
            case 2
           M{i,3}=My{2}(find(sum(abs(repmat(M{i,3},size(y2,1),1)-y2),2)==0));  
            case 3
           M{i,3}=My{3}(find(sum(abs(repmat(M{i,3},size(y3,1),1)-y3),2)==0));
            case 4
           M{i,3}=My{4}(find(sum(abs(repmat(M{i,3},size(y4,1),1)-y4),2)==0));
       end
    end
MD(dd)=0
    for i=1:1:size(M,1)
    MD(dd)=MD(dd)+M{i,1}*M{i,2}*M{i,3};    
    end
    dd=dd+1;
lines=lines+1;
end

fclose(fid);

MDD=cell(size(Mx));
if d==2
   switch deg
    case 1
    MDD{1}=MD(1:2);    
    case 2
    MDD{1}=MD(1:2);
    MDD{2}=MD(3:5);
    case 3
    MDD{1}=MD(1:2);
    MDD{2}=MD(3:5);
    MDD{3}=MD(6:9);
    case 4
    MDD{1}=MD(1:2);
    MDD{2}=MD(3:5);
    MDD{3}=MD(6:9);
    MDD{4}=MD(10:14);
   end
end
 if d==3
   switch deg
    case 1
    MDD{1}=MD(1:3);    
    case 2
    MDD{1}=MD(1:3);
    MDD{2}=MD(4:9);
    case 3
    MDD{1}=MD(1:3);
    MDD{2}=MD(4:9);
    MDD{3}=MD(10:19);
    case 4
    MDD{1}=MD(1:3);
    MDD{2}=MD(4:9);
    MDD{3}=MD(10:19);
    MDD{4}=MD(20:34);
   end
end   
   
end