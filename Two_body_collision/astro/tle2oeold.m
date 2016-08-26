% [oe,epoch,M,E,satname] = tle2oe(fname)
%      fname is a filename string for a file containing
%            a two-line element set (TLE)
%      oe is a 1x6 matrix containing the orbital elements
%            [a e i Om om nu]
% This function returns the orbital elements from TLE.txt
% Calls Newton iteration function file EofMe.m
% SEE ALSO rv2oe, oe2rv
function [oe,epoch,M,E,satname] = tle2oe(fname);

% Open the file up and scan in the elements
fid = fopen(fname, 'rb');
A = fscanf(fid,'%13c%*s',1);
B = fscanf(fid,'%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
C = fscanf(fid,'%d%6d%f%f%f%f%f%f',[1,8]);
fclose(fid);
satname=A;

% The value of mu is for the earth
mu = 3.986012e5;

%Calculate epoch in julian days
epoch = B(1,4);
% ndot = B(1,5);
% n2dot = B(1,6);

% Assign variables to the orbital elements
i = C(1,3)*pi/180;         %Inclination
Om = C(1,4)*pi/180;        %Right Ascension of the Ascending Node
e = C(1,5)/1e7;            %Eccentricity
om = C(1,6)*pi/180;        %Argument of periapsis
M = C(1,7)*pi/180;         %Mean anomaly
n = C(1,8)*2*pi/(24*3600); %Mean motion

% Calculate the semi-major axis
a = (mu/n^2)^(1/3);

% Calculate the eccentric anomaly using Mean anomaly
E = EofMe(M,e,1e-10);
cosnu = (e-cos(E)) / (e*cos(E)-1);
sinnu=((a*sqrt(1-e*e))/(a*(1-e*cos(E))))*sin(E);
nu=atan2(sinnu,cosnu);
if (nu<0), nu=nu+2*pi; end

% Return the orbital elements in a 1X6 matrix
oe = [a e i Om om nu];