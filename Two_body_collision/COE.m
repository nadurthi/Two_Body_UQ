function [sma,ecc,inc,RAAN,aop,nu,TP,hvec,evec] = COE(reci,veci)
global J2 mu Re

hvec = cross(reci, veci);
hmag = norm(hvec);

nvec = cross([0;0;1],hvec);
nmag = norm(nvec);

evec = (((norm(veci))^2-mu/norm(reci)).*reci-dot(reci,veci).*veci)/mu;
ecc = norm(evec);

nrg = (((norm(veci))^2)/2)-mu/norm(reci);

sma = -mu/(2*nrg);

TP = 2*pi*sqrt((sma^3)/mu);

inc = acos(hvec(3)/hmag);

if nvec(2) < 0
    RAAN = 2*pi - acos(nvec(1)/norm(nvec));
else
    RAAN = acos(nvec(1)/norm(nvec));
end

if evec(3) < 0
    aop= 2*pi - acos((dot(nvec,evec))/(ecc*norm(nvec)));
else
    aop =  acos((dot(nvec,evec))/(ecc*norm(nvec)));
end
if dot(reci,veci) < 0
    nu = 2*pi - acos((dot(evec,reci)/(ecc*norm(reci))));
else
    nu = acos((dot(evec,reci)/(ecc*norm(reci))));
end
