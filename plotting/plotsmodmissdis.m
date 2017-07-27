TTstep{1}='14-07-10';
TTstep{2}='14-10-52';
TTstep{3}='14-14-21';
TTstep{4}='14-17-54';
TTstep{5}='14-21-40';
TTstep{6}='14-25-09';
TTstep{7}='14-28-49';
TTstep{8}='14-32-46';
TTstep{9}='14-36-17';
TTstep{10}='14-40-16';
TTstep{11}='14-43-56';
TTstep{12}='14-47-45';
TTstep{13}='14-51-56';
TTstep{14}='14-55-42';
TTstep{15}='14-59-46';
TTstep{16}='15-03-50';

for i=1:1:16
hgload(strcat('MISSxy1_',TTstep{i}))
saveas(gcf,strcat('Missxy1',num2str(49+i)),'pdf')

hgload(strcat('MISSxy2_',TTstep{i}))
view(153,58)
saveas(gcf,strcat('Missxy2',num2str(49+i)),'pdf')


hgload(strcat('MISSxz1_',TTstep{i}))
saveas(gcf,strcat('Missxz1',num2str(49+i)),'pdf')

hgload(strcat('MISSxz2_',TTstep{i}))
view(153,58)
saveas(gcf,strcat('Missxz2',num2str(49+i)),'jpg')
end


%% this is for the paper production of miss- distance pdf over the square f side 2R at the origin

% for i=59:62
i=62
 close all   
 plot1=openfig(strcat('MISSGaussxz1_',num2str(i),'.fig'))
 BB=[get(gca,'XLim') get(gca,'YLim')];
 
 h1 = findall(plot1,'Type','hggroup');
    X1 = get(h1(1),'XData');
    Y1 = get(h1(1),'YData');
    Z1 = get(h1(1),'ZData');
    
    figure
 axesObjs = get(plot1, 'Children');
    dataObjs = get(axesObjs, 'Children');
    objTypes = get(dataObjs, 'Type');
    xdata = get(dataObjs, 'XData');
    ydata = get(dataObjs, 'YData');
    plot(xdata{2},ydata{2},'k','linewidth',2)
%     text(10,2,'\Omega','FontSize',24)
% [1,1]  to [0,0]
x=[0.8,0.71];
y=[0.8,0.75];
    txtar = annotation('textarrow',x,y,'String','\Omega','FontSize',20);
               
    hold on
    
    ls=logspace(-3,-4.5,4);
    [C,h] =contour(X1,Y1,Z1,ls,'linewidth',2);
    
    for kk=1:1:length(ls)
        ind=find(C(1,:)==ls(kk));
        for j=1:1:length(ind)
            %    plot(C(1,ind(j)+1:ind(j)+C(2,ind(j))),C(2,ind(j)+1:ind(j)+C(2,ind(j))))
            XX=C(1,ind(j)+1:ind(j)+C(2,ind(j)));
            YY=C(2,ind(j)+1:ind(j)+C(2,ind(j)));
            N=length(XX);
            hi=1+N*rand;
            hi=round(hi);
            text(XX(hi(1)), YY(hi(1)), num2str(kk),'BackgroundColor',[1 1 1],'FontSize',25)
        end
    end
    
         plot_prop_paper
     grid
     axis([-300,100,-4,1])
% end



