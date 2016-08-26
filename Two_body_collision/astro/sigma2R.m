% R = sigma2R(sigma)
%     sigma is the modified rodrigues parameter vector
%     returns the rotation matrix
%     ref:  HS & JLJ, JAS, 44(1):1, 1996
function R = sigma2R(sigma)
[m,n]=size(sigma);
if (m==1 & n==3)
    sigma=sigma';
end
if ((m~=1 & n~=1) | (m~=3 & n~=3))
   disp('sigma is not a 3 x 1 matrix');
   R=NaN*ones(3,3);
   return;
end
S = 1-sigma'*sigma;
T = (1+sigma'*sigma)^2;
crs = cr(sigma);
R = eye(3)-4*S*crs/T+8*crs*crs/T;