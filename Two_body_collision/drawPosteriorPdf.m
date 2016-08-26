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
%fig=figure;
index1 = 1;
index2 = 2;
sigBound1 = 1;%500 for range
sigBound2 = 10;

rSite = [0;0;0];
indexNumber = [1 200 401 498 699];
indexSize = length(indexNumber);
muVal{1} = num2cell(imu);
sigVal{1} = num2cell(isig);
% limitVal = [-6e-5 6e-5 -6e-5 6e-5;      %for angles-angles
%              -0.04 -0.015 3.25e-3 7e-3;
%              -0.015 0.03 -2e-3 4e-3;
%              0.05 0.12 -3.138 -3.122;
%              0 0.1 -3.138 -3.124];
limitVal = [-6e-5 6e-5;
            -0.035 -0.016;
            -0.015 0.03;
            0.05 0.11;
            0 0.09];

for counter = 3 %last was 301 - 400
    ctr = indexNumber(counter);
    meanTemp = cell2mat(muVal{ctr});
    sigTemp = cell2mat(sigVal{ctr});
%     maxx1 = meanVal(index1) + 5*sqrt(max(max(sigTemp)));
%     minx1 = meanVal(index1) - 5*sqrt(max(max(sigTemp)));
%     maxx2 = meanVal(index2) + 5*sqrt(max(max(sigTemp)));
%     minx2 = meanVal(index2) - 5*sqrt(max(max(sigTemp)));
    cd 'E:\MS@Buffalo\Research\LEO'
    [meanSpherical(:,1), sigmaSpherical] = meanSigmaFromCartesian2Spherical(meanTemp,sigTemp);
    cd 'E:\MS@Buffalo\Research\LEODragEKF'
    maxx1 = meanSpherical(index1,1) + sigBound1*sqrt(max(max(sigmaSpherical)));
    minx1 = meanSpherical(index1,1) - sigBound1*sqrt(max(max(sigmaSpherical)));
    maxx2 = meanSpherical(index2,1) + sigBound2*sqrt(max(max(sigmaSpherical)));
    minx2 = meanSpherical(index2,1) - sigBound2*sqrt(max(max(sigmaSpherical)));
    y = 0;
    limitTemp = limitVal(counter,:);
    
    ctr
    iter = ctr;
    mu = zeros(2,1);
    
    %x1 = linspace(minx1,maxx1,5000);
    x1 = linspace(6985,7020,5000);
    %x2 = linspace(-0.03,0.04,5000);
    x2 = linspace(limitTemp(1),limitTemp(2),5000);
    %x2 = linspace(limitTemp(3),limitTemp(4),5000);
    
        
    [X1,X2] = meshgrid(x1,x2);
    mu(1) = meanSpherical(index1);
    mu(2) = meanSpherical(index2);
    sigTemp = [sigmaSpherical];
    sigTemp = [sigTemp(index1,index1) sigTemp(index1,index2);sigTemp(index2,index1) sigTemp(index2,index2)];
    %sigTemp = cell2mat(sigVal{ctr});
    %sigTemp = [sigTemp(index1,index1) sigTemp(index1,index2);sigTemp(index2,index1) sigTemp(index2,index2)];
    y = mvnpdf([X1(:) X2(:)],[mu(1) mu(2)],sigTemp);
 
    y = reshape(y,length(x2),length(x1));
    
    hold on;
    v = (0.01.*max(max(y)));
   contour(x1,x2,y,[v v],'color',[0.2 0.2 0.6],'LineWidth',2);
   hold on;
   plot(mu(1),mu(2),'*','color',[0.2 0.2 0.6]);
end