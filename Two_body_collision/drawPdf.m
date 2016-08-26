clc
nSteps = 740;
tt=[0:29.173785213265145:6*60*60];
t = tt/60/60;
saveFile3(1,1) = 1;
Dir = 'FPKECKECombinedDrag4Sigma/CKE';

x = numel(time.tspan);

x1 = linspace(-50000,50000,1000);
x2 = linspace(-50000,50000,1000);
        
[X1,X2] = meshgrid(x1,x2);

isig =blkdiag(0.01,0.01,0.01,0.000001,0.000001,0.000001);
imu = [7000 0 0 0 -1.0374090357 7.4771288355];
maxx1 = imu(1) + 1*sqrt(isig(1,1));
minx1 = imu(1) - 1*sqrt(isig(1,1));
maxx2 = imu(2) + 1*sqrt(isig(2,2));
minx2 = imu(2) - 1*sqrt(isig(2,2));

iter = 1;
clc
X = zeros(1,125000);
fig=figure;
index1 = 1;
index2 = 2;
sigBound1 = 1;%500 for range
sigBound2 = 10;
for ctr = 699 %last was 301 - 400
    meanVal = cell2mat(muVal{ctr});
    sigTemp = cell2mat(sigVal{ctr});
    maxx1 = meanVal(index1) + 5*sqrt(max(max(sigTemp)));
    minx1 = meanVal(index1) - 5*sqrt(max(max(sigTemp)));
    maxx2 = meanVal(index2) + 5*sqrt(max(max(sigTemp)));
    minx2 = meanVal(index2) - 5*sqrt(max(max(sigTemp)));
    ctr
    iter = ctr;
    mu = zeros(2,1);
    y = 0;
    
    x1 = linspace(minx1,maxx1,3500);
    x2 = linspace(minx2,maxx2,3500);
        
    [X1,X2] = meshgrid(x1,x2);
    mu(1) = meanVal(index1);
    mu(2) = meanVal(index2);
    sigTemp = cell2mat(sigVal{ctr});
    sigTemp = [sigTemp(index1,index1) sigTemp(index1,index2);sigTemp(index2,index1) sigTemp(index2,index2)];
    y = mvnpdf([X1(:) X2(:)],[mu(1) mu(2)],sigTemp);
 
    y = reshape(y,length(x2),length(x1));
    
    hold on;
    v = (0.01.*max(max(y)));
   contour(x1,x2,y,[v v],'color',[0.2 0.2 0.6],'LineWidth',2);
   hold on;
   plot(mu(1),mu(2),'*','color',[0.2 0.2 0.6]);
   hold on;
end