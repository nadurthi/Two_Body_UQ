plot(XCM(:,2),XCM(:,3),'bo')
XY=[XCM(:,2),XCM(:,3)];
bndl_xy=[min(XY(:,1))-500,min(XY(:,2))-500];
bndu_xy=[max(XY(:,1))+500,max(XY(:,2))+500];
[XYt,dxy]=transform_domain(XY,bndl_xy,bndu_xy,-1*ones(1,2),1*ones(1,2));

plot(XYt(:,1),XYt(:,2),'bo')
axis([-1,1,-1,1])
w=ones(N,1)/N;
 [y1,M1xy]=Cal_moments_samples(XYt,w,1,'central');
 [y2,M2xy]=Cal_moments_samples(XYt,w,2,'central');
 [y3,M3xy]=Cal_moments_samples(XYt,w,3,'central');
 [y4,M4xy]=Cal_moments_samples(XYt,w,4,'central');
 
 mu =M1xy'; 
SIGMA = [M2xy(1) M2xy(2);M2xy(2) M2xy(3)]; 
[x,y]=meshgrid(-1:0.01:1);
p=zeros(size(x));
for i=1:1:size(x,1)
    for j=1:1:size(x,2)
p(i,j) = mvnpdf([x(i,j) y(i,j)],mu,SIGMA); 
    end
end
contour(x,y,p,[linspace(0,0.1,10),linspace(1,20,10)]) 
hold on
plot(XYt(:,1),XYt(:,2),'bo')
xlabel('y')
ylabel('z')
axis([-1,1,-1,1])
plot_prop_paper


XY=[XIR(:,2),XIR(:,3)];
% bndl_xy=[min(XY(:,1))-500,min(XY(:,2))-500];
% bndu_xy=[max(XY(:,1))+500,max(XY(:,2))+500];
[XYt,dxy]=transform_domain(XY,bndl_xy,bndu_xy,-1*ones(1,2),1*ones(1,2));

% plot(XYt(:,1),XYt(:,2),'bo')
% axis([-1,1,-1,1])
w=ones(N,1)/N;
 [y1,M1xy]=Cal_moments_samples(XYt,w,1,'central');
 [y2,M2xy]=Cal_moments_samples(XYt,w,2,'central');

 
 mu =M1xy'; 
SIGMA = [M2xy(1) M2xy(2);M2xy(2) M2xy(3)]; 
[x,y]=meshgrid(-1:0.01:1);
p2=zeros(size(x));
for i=1:1:size(x,1)
    for j=1:1:size(x,2)
p2(i,j) = mvnpdf([x(i,j) y(i,j)],mu,SIGMA); 
    end
end
contour(x,y,p2,[linspace(0,0.1,10),linspace(1,20,10)]) 
hold on
plot(XYt(:,1),XYt(:,2),'ro')
axis([-1,1,-1,1])