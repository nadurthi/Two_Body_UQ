% [lat,long] = oe2ssp(oe,dt,gst,mu)
%
%  Author:  Capt Ralph Sandfry                  23 Feb 99
%
%  This function calculates the latitude and longitude of
%  a satellite's sub-satellite point at epoch and at other  
%  specified times following epoch.  The only mandetory
%  input is the orbital element vector,oe.  If dt is omitted
%  then the lat/long at epoch is returned.  If gst is omitted
%  then the local sidereal time is returned instead of
%  longitude.  Earth's gravitational parameter is the default.
%
%  INPUTS:
%       oe      - vector of the 6 orbital elements,
%                 a (km),e,i,Om,om,nu (radians) 
%       dt      - vector of time increments from epoch    sec
%       gst0    - Greenwich Sidereal Time at epoch        0.0 to 2Pi rad
%       mu      - Gravitational Parameter                 km^3/s^2
%
%  OUTPUTS:
%       lat     -  latitude of sub-satellite point(s)     rad
%       long    -  longitude of sub-satellite point(s)    rad
%                  Note:  if gst=0, 
%                            long = local sidereal time
%  LOCALS:
%       num     -  length of time vector
%       I       -  counter
%       rv      -  position vector                         km
%       vv      -  velocity vector                         km/s
%       oe_dt   -  orbital elements for time dt
%       gst     -  Greenwich Sidereal Time                 rad
%       param   -  gst-local sidereal time switch          0 or 1
%       
%  CONSTANTS:
%       twopi   -  2*pi
%       omega   -  Earth's rotation rate                   rad/s
%
%  COUPLING:
%       oeofdt  -  propogates orbital elements forward in time
%       oe2rv   -  calculates position and velocity vectors from
%                  orbital elements
%  REFERENCES:
%       AOE 4984 Class Notes, Dr C.D. Hall 
%       Spring 1999,Virginia Tech 
%
%
function [lat,long] = oe2ssp(oe,dt,gst0,mu)

if nargin < 4  
   mu = 398600.5;     
end

param = 1;

if nargin < 3
   gst0 = 0;     %  function will return local sidereal time
   param = 0;    %  in place of longitude
end              

if nargin < 2
   dt = 0;
end

num = length(dt);
omega = 7.2921159e-5;

for I = 1:num
   
   oe_dt = oeofdt(oe,dt(I),mu);            % update orbital elements
   
   [rv,vv] = oe2rv(oe_dt,mu);              % calculate position & velocity from oe
   
   gst(I) = gst0 + param*omega*dt(I);      % update GST from epoch
   
   long(I) = atan2(rv(2),rv(1))-gst(I);    % sub-satellite point, longitude
   while long(I)<0
      long(I)=long(I)+2*pi;
   end
      
   lat(I) =  asin(rv(3)/norm(rv));         % sub-satellite point, latitude
end


long = long';
lat = lat';