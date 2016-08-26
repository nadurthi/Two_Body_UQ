% oe = tle2oee(fname)
%      fname is a filename string for a file containing
%            a two-line element set (TLE)
%      oe is a 1x6 matrix containing the orbital elements
%            [a e i Om om E]
% This function returns the orbital elements from TLE.txt
% Calls Newton iteration function file EofMe.m

function oe = tle2oe()

% Open the file up and scan in the elements
fid = fopen('egtles.m', 'rb');
A = fscanf(fid,'%11c%*s',1);
B = fscanf(fid,'%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
C = fscanf(fid,'%d%6d%f%f%f%f%f%f',[1,8]);
fclose(fid);

% The value of mu is for the earth
mu = 3.98601e5;

%Calculate epoch
epoch = B(1,4)*24*3600;

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
% Call the function newton2.m
E = EofMe(M,e,1e-10);

% Return the orbital elements in a 1X6 matrix
oe = [a e i Om om E];
