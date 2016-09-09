r1 = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
opt = odeset('reltol',1e-12,'abstol',1e-12);

N=5000;
 T=0:30:172800;
P0 =blkdiag(0.01,0.01,0.01,1e-009,1e-009,1e-009);
X1mc0=mvnrnd(r1,P0,N);
X1mc=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,X1mc0(i,:)',opt);
   X1mc(:,:,i)=X;
end