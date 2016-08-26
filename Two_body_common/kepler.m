function [E]=kepler(M, e, E0, tol)
% Code to solve the Kepler's Equation
% M - (E - e*sin(E)) = 0 
% M - Mean anomaly
% e - orbit eccentricity
% E0 - Eccentric anomaly (initial guess)
% E - Eccentric anomaly
% tol - tolerance for Newton's Root Solver

nitr = 10; E = E0;

for k = 1:nitr
    fx = M - (E-e*sin(E));
    dfx = -(1-e*cos(E));
    
    dE = -fx/dfx;
    
    if (abs(dE)>tol)
         k;
         E = E + dE;
    else
        disp(['Solution Converged in ',num2str(k),' iterations']);
        break;
    end
end
    