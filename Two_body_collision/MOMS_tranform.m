function MDD=MOMS_tranform(Mx,A)
d=length(Mx{1});
deg=size(Mx,1);
MDD=cell(size(Mx));
[X,w]=GH_points(zeros(d,1),eye(d),4);
[y1,M]=Cal_moments_samples(X,w,1,'raw');
[y2,M]=Cal_moments_samples(X,w,2,'raw');
[y3,M]=Cal_moments_samples(X,w,3,'raw');
[y4,M]=Cal_moments_samples(X,w,4,'raw');

if d==3
fid = fopen('Moms_transform_txt.m');
totdeg=25;
MDD{1}=(A*Mx{1}')';
ss=A*[Mx{2}(1),Mx{2}(3),Mx{2}(2);Mx{2}(3),Mx{2}(6),Mx{2}(4);Mx{2}(2),Mx{2}(4),Mx{2}(5)]*A';
MDD{2}=[ss(1,1),ss(1,3),ss(1,2),ss(2,3),ss(3,3),ss(2,2)];
end
if d==2
fid = fopen('Moms_transform_txt2D.m');    
totdeg=9;
MDD{1}=(A*Mx{1}')';
ss=A*[Mx{2}(1),Mx{2}(2);Mx{2}(2),Mx{2}(3)]*A';
MDD{2}=[ss(1,1),ss(1,2),ss(2,2)];
end

lines=1;
dd=1;
MD=0;
while lines<=totdeg
    tline = fgets(fid);
    tline=strcat(tline,',');
    ind=strfind(tline,',');
    M=cell(length(ind)-1,3);
    for i=1:1:length(ind)-1
        ss=tline(ind(i)+1:ind(i+1)-1);
        ss=strcat(ss,',');
        M{i,1}=getNUMcoeff(ss);
        if d==2
        M{i,2}=prod(prod(A.^[getSTRpower(ss,'a11'),getSTRpower(ss,'a12');getSTRpower(ss,'a21'),getSTRpower(ss,'a22')]));
        M{i,3}=[getSTRpower(ss,'x1'),getSTRpower(ss,'x2')];    
        else
        M{i,2}=prod(prod(A.^[getSTRpower(ss,'a11'),getSTRpower(ss,'a12'),getSTRpower(ss,'a13');getSTRpower(ss,'a21'),getSTRpower(ss,'a22'),getSTRpower(ss,'a23');getSTRpower(ss,'a31'),getSTRpower(ss,'a32'),getSTRpower(ss,'a33')]));
        M{i,3}=[getSTRpower(ss,'x1'),getSTRpower(ss,'x2'),getSTRpower(ss,'x3')];
        end
        switch sum(M{i,3})
           case 1
           M{i,3}=Mx{1}(find(sum(abs(repmat(M{i,3},size(y1,1),1)-y1),2)==0));    
            case 2
           M{i,3}=Mx{2}(find(sum(abs(repmat(M{i,3},size(y2,1),1)-y2),2)==0));  
            case 3
           M{i,3}=Mx{3}(find(sum(abs(repmat(M{i,3},size(y3,1),1)-y3),2)==0));
            case 4
           M{i,3}=Mx{4}(find(sum(abs(repmat(M{i,3},size(y4,1),1)-y4),2)==0));
       end
    end
    MD(dd)=0;
    for i=1:1:size(M,1)
        MD(dd)=MD(dd)+M{i,1}*M{i,2}*M{i,3};      
    end
    dd=dd+1;
    lines=lines+1;
end

fclose(fid);

if d==2
   switch deg
    
    case 3
    MDD{3}=MD(1:4);
    case 4
    MDD{3}=MD(1:4);
    MDD{4}=MD(5:9);
   end
end
 if d==3
   switch deg
    
    case 3
    MDD{3}=MD(1:10);
    case 4
    MDD{3}=MD(1:10);
    MDD{4}=MD(11:25);
   end
end   
end


function p=getSTRpower(str,var)

ind=strfind(str,var);
if length(ind)==0
    p=0;
    return;
end

if strcmp(str(ind+length(var)),'^')==1
    p=str2num(str(ind+length(var)+1));
else
    p=1;
end


end


function n=getNUMcoeff(str)
str(strfind(str,' '))=[];
ind=strfind(str,'a');
str(ind(1):end)=[];
if length(str)==1 || length(str)==0
   n=1;
else
    n=str2num(str);
end
end