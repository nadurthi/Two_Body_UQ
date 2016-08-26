function [T0,TT0,YT,DT,n0d,n0dd,Bstr,EphType,Type,i,e,omg,Omg,M,n]=orb_element(L1,L2)
%% taken from celstrak matlab files
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(72);
deg2rad  =   pi / 180.0;         %  0.01745329251994330;  % [deg/rad]
    xpdotp   =  1440.0 / (2.0*pi);   % 229.1831180523293;  % [rev/day]/[rad/min]   
    revnum = 0; 
    elnum  = 0;
    year   = 0; 
    satrec.error = 0;
    
    for (j = 11:16)
        if (longstr1(j) == ' ')
            longstr1(j) = '_';
        end
    end

    if (longstr1(45) ~= ' ')
        longstr1(44) = longstr1(45);
    end
    longstr1(45) = '.';
     
    if (longstr1(8) == ' ')
        longstr1(8) = 'U';
    end

    if (longstr1(10) == ' ')
        longstr1(10) = '.';
    end

    for (j = 46:50)
        if (longstr1(j) == ' ')
            longstr1(j) = '0';
        end
    end
    if (longstr1(52) == ' ')
        longstr1(52) = '0';
    end
    if (longstr1(54) ~= ' ')
        longstr1(53) = longstr1(54);
    end
    longstr1(54) = '.';

    longstr2(26) = '.';
     
    for (j = 27:33)
        if (longstr2(j) == ' ')
            longstr2(j) = '0';
        end
    end
     
    if (longstr1(63) == ' ')
        longstr1(63) = '0';
    end

    if ((length(longstr1) < 68) || (longstr1(68) == ' '))
        longstr1(68) = '0';
    end

    % parse first line
    carnumb = str2num(longstr1(1));
    satrec.satnum = str2num(longstr1(3:7));
    classification = longstr1(8);
    intldesg = longstr1(10:17);
    satrec.epochyr = str2num(longstr1(19:20));
    satrec.epochdays = str2num(longstr1(21:32));
    satrec.ndot = str2num(longstr1(34:43));
    satrec.nddot = str2num(longstr1(44:50));
    nexp = str2num(longstr1(51:52));
    satrec.bstar = str2num(longstr1(53:59));
    ibexp = str2num(longstr1(60:61));
    numb = str2num(longstr1(63));
    elnum = str2num(longstr1(65:68));
    cardnumb = str2num(longstr2(1));
    % parse second line
        satrec.satnum = str2num(longstr2(3:7));
        satrec.inclo = str2num(longstr2(8:16));
        satrec.nodeo = str2num(longstr2(17:25));
        satrec.ecco = str2num(longstr2(26:33));
        satrec.argpo = str2num(longstr2(34:42));
        satrec.mo = str2num(longstr2(43:51));
        satrec.no = str2num(longstr2(52:63));
        revnum = str2num(longstr2(64:68));
    
        
        satrec.no   = satrec.no / xpdotp; %//* rad/min
    satrec.nddot= satrec.nddot * 10.0^nexp;
    satrec.bstar= satrec.bstar * 10.0^ibexp
    
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     *deg2rad;
    
    if (satrec.epochyr < 57)
         year= satrec.epochyr + 2000;
       else
         year= satrec.epochyr + 1900;
     end;
     
     
     [mon,day,hr,minute,sec] = days2mdh ( year,satrec.epochdays );
     satrec.jdsatepoch = jday( year,mon,day,hr,minute,sec );
     
     %     // ---- convert to sgp4 units ----
    satrec.a    = (satrec.no*tumin)^(-2/3);                % [er]
    satrec.ndot = satrec.ndot  / (xpdotp*1440.0);          % [rad/min^2]
    satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);     % [rad/min^3]
 %%   
  


YT=str2double(L1(19:20));
DT=str2double(L1(21:32));
n0d=str2double(L1(34:43));
n0dd=str2double(L1(45:52));
Bstr=str2double(L1(54:61));
EphType=str2double(L1(63));
% 
% i=str2double(L2(9:16))*pi/180;
% Omg=str2double(L2(18:25))*pi/180; %Right Ascension of ascending node/ lon
% omg=str2double(L2(35:42))*pi/180;
% e=str2double(L2(27:33));
% M=str2double(L2(44:51))*pi/180;%mean anamoly
% n=str2double(L2(53:63));%mean motion revs/day
% n=n*(2*pi)/(24*60*60); %radians/s;
% 

if EphType<=1
    Type='SGP';
end
if EphType==2
    Type='SGP4';
end
if EphType==3
    Type='SDP4';
end
if EphType==4
    Type='SGP8';
end
if EphType==5
    Type='SDP8';
end
% i0=str2double(L2(9:16));% degrees
% W0=str2double(L2(18:25)); %right ascension of the scending node degree: ascending node at epoch
% e0=str2double(L2(27:33));
% w0=str2double(L2(35:42));% arg at perigee degrees
% M0=str2double(L2(44:51));
% n0=str2double(L2(53:63)); %mean motion revs per day


i=str2double(L2(9:16))*pi/180;
Omg=str2double(L2(18:25))*pi/180; %Right Ascension of ascending node/ lon
omg=str2double(L2(35:42))*pi/180;
e=str2double(L2(27:33));
M=str2double(L2(44:51))*pi/180;%mean anamoly
n=str2double(L2(53:63));%mean motion revs/day
n=n*(2*pi)/(24*60*60); %radians/s;

if YT>=0 && YT<=56
Y=YT+2000;
end
if YT>=57 && YT<=99
Y=YT+1900;
end
D=floor(DT);
h=(DT-D)*24;
m=(h-floor(h))*60;
s=(m-floor(m))*60;
T0=[Y,D,h,m,s];
vec_day=zeros(1,12);
for i=1:1:12
    if i==1
        vec_day(i)=eomday(Y, i);
    else
    
vec_day(i)=vec_day(i-1)+eomday(Y, i);
    end
end
ind=find(vec_day>=D);
month=ind(1);
TT0=[floor(Y),month,D-vec_day(month-1),h,m,s];



