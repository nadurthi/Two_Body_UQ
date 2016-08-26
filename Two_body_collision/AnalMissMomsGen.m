Y1=[     1     0     0
0     1     0
0     0     1];
%
Y2=  [2     0     0
1     0     1
1     1     0
0     1     1
0     0     2
0     2     0];
%
Y3=[ 3     0     0
2     0     1
2     1     0
1     1     1
1     0     2
1     2     0
0     2     1
0     1     2
0     0     3
0     3     0];
%
Y4=[ 4     0     0
3     0     1
3     1     0
2     1     1
2     0     2
2     2     0
1     2     1
1     1     2
1     0     3
1     3     0
0     3     1
0     2     2
0     1     3
0     0     4
0     4     0];
%
Y5=[ 5     0     0
4     0     1
4     1     0
3     1     1
3     0     2
3     2     0
2     2     1
2     1     2
2     0     3
2     3     0
1     3     1
1     2     2
1     1     3
1     0     4
1     4     0
0     4     1
0     3     2
0     2     3
0     1     4
0     0     5
0     5     0];
%
Y6=[6     0     0
5     0     1
5     1     0
4     1     1
4     0     2
4     2     0
3     2     1
3     1     2
3     0     3
3     3     0
2     3     1
2     2     2
2     1     3
2     0     4
2     4     0
1     4     1
1     3     2
1     2     3
1     1     4
1     0     5
1     5     0
0     5     1
0     4     2
0     3     3
0     2     4
0     1     5
0     0     6
0     6     0];


%%
syms x1 x2 x3 y1 y2 y3
d1=x1-y1;
d2=x2-y2;
d3=x3-y3;
fid = fopen('AnalMissMoms.m', 'w');
for i=1:1:size(Y1,1)
fprintf(fid,'%s ;\n',strcat('M1(',num2str(i),')=',char(expand(d1^Y1(i,1)*d2^Y1(i,2)*d3^Y1(i,3)))))
end
for i=1:1:size(Y2,1)
fprintf(fid,'%s ;\n',strcat('M2(',num2str(i),')=',char(expand(d1^Y2(i,1)*d2^Y2(i,2)*d3^Y2(i,3)))))
end
for i=1:1:size(Y3,1)
fprintf(fid,'%s ;\n',strcat('M3(',num2str(i),')=',char(expand(d1^Y3(i,1)*d2^Y3(i,2)*d3^Y3(i,3)))))
end
for i=1:1:size(Y4,1)
fprintf(fid,'%s ;\n',strcat('M4(',num2str(i),')=',char(expand(d1^Y4(i,1)*d2^Y4(i,2)*d3^Y4(i,3)))))
end
for i=1:1:size(Y5,1)
fprintf(fid,'%s ;\n',strcat('M5(',num2str(i),')=',char(expand(d1^Y5(i,1)*d2^Y5(i,2)*d3^Y5(i,3)))))
end
for i=1:1:size(Y6,1)
fprintf(fid,'%s ;\n',strcat('M6(',num2str(i),')=',char(expand((d1^Y6(i,1))*(d2^Y6(i,2))*(d3^Y6(i,3))))))
end
fclose(fid)

%%
fid = fopen('AnalMissMoms6D.m', 'r');
fid2 = fopen('AnalMissMoms6D.m_converted.m', 'w');
lne=fgets(fid);
slin=1;
while(length(lne)>0)
    str=lne;
    i=1;
    ee=zeros(1,3);
    while i<length(lne)
        if lne(i)=='x'
            if lne(i+2)=='^'
                ee(str2num(lne(i+1)))=str2num(lne(i+3));
                i=i+4;
            else
                ee(str2num(lne(i+1)))=1;
                i=i+2;
            end
            %             ee
        else
            i=i+1;
        end
        if (lne(i)=='+' || lne(i)==';' || lne(i)=='-') && sum(ee)>0
            tend=i-1;
            j=tend;
            while lne(j)~='+' && lne(j)~='-' && lne(j)~='='
               j=j-1; 
            end
            tstart=j+1;
            one_lne=lne(1:j);
            midlne=lne(tstart:tend);
            thrd_lne=lne(tend+1:end);
            ind=strfind(midlne,'x');
            x1=ind(1);
            lastx=ind(end);
            j=lastx;
            while midlne(j)~='y' && midlne(j)~=' ' && midlne(j)~='*'
                j=j+1;
            end
            x2=j-1;
             MM=MomentVecorder(ee);
             midlnep=midlne(1:x1-1);
             midlnen=midlne(x2+1:end);
             midlne=strcat(midlnep,strcat('MX',num2str(MM(1)),'(',num2str(MM(2)),')'),midlnen); 
%             midlne = regexprep(midlne, midlne(x1:x2), strcat('MM',num2str(MM(1)),'(',num2str(MM(2)),')'));
            lne=strcat(one_lne,midlne,thrd_lne);
            ee=zeros(1,3);
            i=1;
        end
    end
    fprintf(fid2,'%s \n',lne);
    slin=slin+1
    if slin> 925
        break;
    end
    
    %     keyboard
    %
    lne=fgets(fid);
end
fclose(fid)
fclose(fid2)


%% parsing y
fid = fopen('AnalMissMoms6D.m_converted.m', 'r');
fid2 = fopen('AnalMissMoms6D.m_converted2.m', 'w');
lne=fgets(fid);
slin=1;
while(length(lne)>0)
    str=lne;
    i=1;
    ee=zeros(1,3);
    while i<length(lne)
        if lne(i)=='y'
            if lne(i+2)=='^'
                ee(str2num(lne(i+1)))=str2num(lne(i+3));
                i=i+4;
            else
                ee(str2num(lne(i+1)))=1;
                i=i+2;
            end
            %             ee
        else
            i=i+1;
        end
        if (lne(i)=='+' || lne(i)==';' || lne(i)=='-') && sum(ee)>0
            tend=i-1;
            j=tend;
            while lne(j)~='+' && lne(j)~='-' && lne(j)~='='
               j=j-1; 
            end
            tstart=j+1;
            one_lne=lne(1:j);
            midlne=lne(tstart:tend);
            thrd_lne=lne(tend+1:end);
            ind=strfind(midlne,'y');
            y1=ind(1);
%             lasty=ind(end);
%             j=lasty;
%             while midlne(j)~='y' && midlne(j)~=' ' && midlne(j)~='*'
%                 j=j+1;
%             end
%             x2=j-1;
             MM=MomentVecorder(ee);
             midlnep=midlne(1:y1-1);
%              midlnen=midlne(x2+1:end);
             midlne=strcat(midlnep,strcat('MY',num2str(MM(1)),'(',num2str(MM(2)),')')); 
%             midlne = regexprep(midlne, midlne(x1:x2), strcat('MM',num2str(MM(1)),'(',num2str(MM(2)),')'));
            lne=strcat(one_lne,midlne,thrd_lne);
            ee=zeros(1,3);
            i=1;
        end
    end
    fprintf(fid2,'%s \n',lne);
    slin=slin+1
    if slin> 925
        break;
    end
    
    %     keyboard
    %
    lne=fgets(fid);
end
fclose(fid)
fclose(fid2)