clc
nSteps = 740;
tt=[0:29.173785213265145:6*60*60];
t = tt/60/60;

x = numel(time.tspan);
prior.n = 49;
x1 = linspace(-50000,50000,1000);
x2 = linspace(-50000,50000,1000);
        
[X1,X2] = meshgrid(x1,x2);

isig =blkdiag(0.01,0.01,0.01,0.000001,0.000001,0.000001);
imu = [7000 0 0 0 -1.0374090357 7.4771288355];
maxx1 = imu(1) + 1*sqrt(isig(1,1));
minx1 = imu(1) - 1*sqrt(isig(1,1));
maxx2 = imu(2) + 1*sqrt(isig(2,2));
minx2 = imu(2) - 1*sqrt(isig(2,2));
Dir = 'plotForMeasurementUpdate/measUpdateAnglesOnly';
iter = 1;
clc
%plot(pdfVal(1),pdfVal(2),'*','color',[0.5 0.5 0.2]);
index1 = 1;
index2 = 2;
sigBound1 = 100;
sigBound2 = 100;
for ctr = 699 %last was 301 - 400
    sigTemp = sigForMeas;
    meanValue = meanForMeas;
    maxx1 = meanValue(index1) + sigBound1*sqrt(max(max(sigTemp)));
    minx1 = meanValue(index1) - sigBound1*sqrt(max(max(sigTemp)));
    maxx2 = meanValue(index2) + sigBound2*sqrt(max(max(sigTemp)));
    minx2 = meanValue(index2) - sigBound2*sqrt(max(max(sigTemp)));
    ctr
    iter = ctr;
    %mu = zeros(prior.n,2);
    y = 0;
    
%     x1 = linspace(minx1,maxx1,5000);
%     x2 = linspace(minx2,maxx2,5000);
    %x1 = linspace(8100,8400,1000);
    x1 = linspace(-0.50,-0.40,1000);
    x2 = linspace(2.44,2.46,1000);
        
    [X1,X2] = meshgrid(x1,x2);
    mu(1) = meanValue(index1);
    mu(2) = meanValue(index2);
    sigTemp = [sigTemp(index1,index1) sigTemp(index1,index2);sigTemp(index1,index2) sigTemp(index2,index2)];
    y = mvnpdf([X1(:) X2(:)],[mu(1) mu(2)],sigTemp);
    y = reshape(y,length(x2),length(x1));
    
    hold on;
    v = (0.01.*max(max(y)));
   %contour(x1,x2,y,[v v],'-g','LineWidth',2);
   contour(x1,x2,y,[v v],'color',[0.7 0.2 0.5],'LineWidth',2);
   hold on;
   plot(mu(1),mu(2),'*','color',[0.7 0.2 0.5]);
   %plot(pdfVal(1),pdfVal(2),'*g');
   set(gca,'fontsize',16);
   xlabel('x(in km)'); ylabel('azimuth(in rad)');
   saveas(gcf,[Dir 'FPKE_XY_' num2str(mcIndex(indexCounter)) '.eps']);
   saveas(gcf,[Dir 'FPKE_XY_' num2str(mcIndex(indexCounter)) '.fig']);
   saveas(gcf,[Dir 'FPKE_XY_' num2str(mcIndex(indexCounter)) '.png']);
   close all
end