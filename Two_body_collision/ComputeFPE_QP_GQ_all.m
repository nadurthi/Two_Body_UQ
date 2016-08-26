function [MM, NN] = ComputeFPE_QP_GQ_all(weig, mu, sig, model, time, k)

% compute matrix MM
MM = zeros(length(weig));
MMtmp = MM;
for i = 1 : length(weig)
    for j = i+1 : length(weig)
        MMtmp(i,j) = getLikelihood(mu{i}-mu{j}, sig{i}+sig{j});
    end
    MM(i,i) = 1/sqrt(det(4.*pi.*sig{i}));
end
MM = MM + MMtmp + MMtmp';  
MM = MM/time.dt^2;

% Generate Collocation Points using Gaussian Quadrature with quadrature
% points for the entire domain
tmp = zeros(model.fn, 2*model.fn*length(weig));
for j = 1 : length(weig)
    for i = 1 : model.fn
        tmpS = chol(sig{j})';
%         tmpBound = myfilter.integration_GQ_bound;
        tmpBound = 3;
        tmp((j-1)*2*model.fn + i) = mu{j} - tmpBound * tmpS(:,i);
        tmp((j-1)*2*model.fn + i+model.fn) = mu{j} + tmpBound * tmpS(:,i);
    end;
end;

xl = min(tmp'); xu = max(tmp');
% npts = myfilter.integration_GQ_no_points;
npts = 1000;
xpts = ones(1,model.fn) * npts;
[xint,wint] = get_colocation(xpts, xl, xu);        

% compute matrix NN
NN = zeros(length(weig));    
for ct = 1:length(wint)
    [Pk, Lk] = ComputeErrCoeff_QP(xint(:,ct), weig, mu, sig, model, time, k);
    NN = NN + wint(ct) * Pk*Lk';
end
