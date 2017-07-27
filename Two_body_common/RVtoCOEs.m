%  function[A] = RVtoCOEs(R,V)
%
%      R = [Ri, Rj, Rk] (radius vector)
%      V = [Vi, Vj, Vk] (velocity vector)
%
%      A = [a, (semi major axis)
%           e, (eccentricity)
%           i, (inclination)
%           node, (RAAN)
%           arg, (argument of perigee)
%           mean (true annomaly) ]
%

function[A] = RVtoCOEs(R,V)

MU=398601.2;

energy=norm(V)^2/2 - MU/norm(R);
a = -MU/2/energy;

H = cross(R,V);
H=H(:)';
E = cross(V,H)/MU - R/norm(R);

inclination = acos(H(3)/norm(H));


   k = [0,0,1];
   N = cross(k,H);
   node = acos(N(1)/norm(N));
   if N(2)<0
      node=2*pi-node;
   end;
   

   arg = acos(  dot(N,E) / ( norm(N)*norm(E) )  );
   if E(3)<0 
      arg = 2*pi - arg;
   end;


   true = acos(dot(E,R)/(norm(E)*norm(R)));
   if dot(R,V)<0
      true = 2*pi - true;
   end;

   
   eccentric = acos( (norm(E)+cos(true))/(1+norm(E)*cos(true)) );
   if true>pi
      eccentric = 2*pi - eccentric;
   end;


   mean = eccentric - norm(E)*sin(eccentric);


   A = [a,
        norm(E),
        (inclination),
        (node),
        (arg),
        (true)];