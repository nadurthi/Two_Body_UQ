
% tsince is in min

% first initialise the orbital elements and find the times to progate to
% and from
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
   global opsmode
   opsmode= 'a';
global idebug dbgfile
endtime.year=2009;
endtime.month=2;
endtime.day=10;
endtime.hr=16;
endtime.min=10;
endtime.sec=0;
dtmin=10;% deltat in min (time step for propagation)

[satrec, startmfe, stopmfe, deltamin] = twoline2rv_modified(72,L1,L2,'m','e',endtime,dtmin);


 %               // call the propagator to get the initial state vector value
            [satrec, ro ,vo] = sgp4 (satrec,  0.0);
            
            
            
 % propagate the shit   
 for tsince = startmfe:deltamin:stopmfe
 [satrec, r, v] = sgp4(satrec,tsince);
 
 end