function [Pk, Lk] = ComputeErrCoeff_QP(x, weig, mu, sig, model, time, k)

%%%
% weig is a N x 1 vector of Gaussian weights
% mu is a Nx1 vector of Gassian Means
% sig is a Nx1 vector of Gaussian Variances
%%%

Pk = zeros(length(weig),1);
Lk = zeros(length(weig),1);
for ct = 1 : length(weig)
    [pg, dpgdx, ddpgdx, dpgdmu, dpgdsig] = GaussianPDF(x, mu{ct}, sig{ct});

    myX = [mu{ct}; reshape(sig{ct}, numel(sig{ct}), 1)];
    myDX = EKF_propagation(time.tspan(k), myX, model);
    mudot = myDX(1:model.fn);
    sigdot = reshape(myDX(model.fn+1:end), model.fn, model.fn);
        
    dpgidt = dpgdmu'*mudot + trace(dpgdsig*sigdot);
    
    fx = feval(model.fx, time.tspan(k), x);
    Jx = feval([model.fx '_jac'], time.tspan(k), x);
    dopgi_dot = - dpgdx' * fx - pg*trace(Jx) + 1/2*trace(model.Q * ddpgdx);
    
    Pk(ct) = pg/time.dt;
    Lk(ct) = dpgidt - dopgi_dot - pg/time.dt;
end;
 
 
 
 
    
