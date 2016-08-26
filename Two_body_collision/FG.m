function [r_tFG, rdot_tFG] = FG(r_0, rdot_0, mu)

r0 = norm(r_0); %% Magnitude of Given Initial Position %%
v0 = norm(rdot_0); %% Magnitude of Given Initial Velocity %%
h = cross(r_0, rdot_0); %% Compute Angular Momentum %%
sigma0 = (r_0'*rdot_0)/sqrt(mu); %% Compute Sigma_0 %%
r_tFG(:,1) = r_0; %% Store Initial Position Vector %%
rdot_tFG(:, 1) = rdot_0; %% Store Initial Velocity Vector %%
t0 = 0; 

for t = 10:10:28800 
    tol=0.001;
    deltaE=1; %% Arbitrary Initial deltaE %%
    Etilde=0;
    while deltaE>tol
        Etildenew = Etilde-(Etilde-(1-r0/a)*sin(Etilde)-(sigma0/sqrt(a))*(cos(Etilde)-1)-sqrt(mu/a^3)*(t-t0))/(1-(1-r0/a)*cos(Etilde)+(sigma0/sqrt(a^3))*sin(Etilde));
        deltaE=abs(Etilde-Etildenew);
        Etilde = Etildenew;
    end
    r = a+(r0-a)*cos(Etildenew)+sigma0*sqrt(a)*sin(Etildenew);
    F = 1-(a/r0)*(1-cos(Etildenew));
    Fdot = (-sqrt(mu*a)/(r*r0))*sin(Etildenew);
    G = (t-t0)+sqrt(a^3/mu)*(sin(Etildenew)-Etildenew);
    Gdot = 1-(a/r)*cos(Etildenew);
    r_tFG(:,t/10+1) = F*r_0+G*rdot_0;
    rdot_tFG(:,t/10+1) = Fdot*r_0+Gdot*rdot_0;
end