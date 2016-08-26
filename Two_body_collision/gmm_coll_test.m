%% GM collision testing
R=0.1;

g.mu=[-2,0;-2,2];
g.P=[reshape([1,0.8;0.8,1],1,4);reshape(0.5*[1,0;0,1],1,4)];
g.w=[0.2;0.8];

h.mu=[2,0;2,2];
h.P=[reshape([1,-0.8;-0.8,1],1,4);reshape(0.5*[1,0;0,1],1,4)];
h.w=[0.2;0.8];

objg=GMM2obj(g);
objh=GMM2obj(h);
N=5000000;
Yg = random(objg,N);
Yh = random(objh,N);
prob_mc=samples_prob_coll(Yg,Yh,R,1)

plot(Yg(1:1e5,1),Yg(1:1e5,2),'ro',Yh(1:1e5,1),Yh(1:1e5,2),'bo')

D=prob_coll_GM(g,h,R)

N2=10000;
probgh=prob_coll_norms(Yg(1:N2,:),ones(N2,1)/N2,Yh(1:N2,:),ones(N2,1)/N2,R)