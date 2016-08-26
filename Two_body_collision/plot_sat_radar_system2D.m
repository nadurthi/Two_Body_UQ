function plot_sat_radar_system2D(XsigSat,k,Radmodel,MeasPairs,yplottruth)
Re=6378.1;



%% plot true sat trajs and current sat positions 
for i=1:1:Radmodel.Nsat
    mu=XsigSat{i,1}(k,1:2);
    
xx1=yplottruth{i,1}(:,1);
xx2=yplottruth{i,1}(:,2);





plot(xx1,xx2,'r:', 'linewidth',0.5)
hold on
if isempty(find(MeasPairs{k}(:,1)==i))
plot(mu(1),mu(2),'ro', 'MarkerFaceColor','r','MarkerSize',6)
else
    plot(mu(1),mu(2),'bo', 'MarkerFaceColor','b','MarkerSize',8)
end
% Create textarrow
% [xaf,yaf] = ds2nfu([mu(1),mu(1)-0.1*mu(1)],[mu(2),mu(2)-0.1*mu(2)]);
% annotation('textarrow',xaf,...
%     yaf,'TextEdgeColor','none',...
%     'String',{'Sat ',num2str(i)});
end
%% plot radar pos and their cones 
for i=1:1:Radmodel.Nrad
   [v,Rot,p]=vec_radar_coordchange2D([0 0],i,Radmodel.RadPos,'local2ecef');
   Rot=Rot';
   p';
plot(v(1),v(2),'ko', 'MarkerFaceColor','k','MarkerSize',6)
if isempty(find(MeasPairs{k}(:,2)==i))
    continue
end

% x1 = ecef2enu([1,0,0],Radmodel.pos(i,:)*1e3);
% x2 = ecef2enu([0,1,0],Radmodel.pos(i,:)*1e3);
% x3 = ecef2enu([0,0,1],Radmodel.pos(i,:)*1e3);
% Rot=[x1(:),x2(:),x3(:)];
% Rot=Rot';
% 
% %% have to shift and rotate the cones
% keyboard
coneang=Radmodel.SensParas(i,1);
R=Radmodel.SensParas(i,2);
% r=linspace(0,R*tan(coneang),30);
% theta=linspace(0,2*pi,30);
% [r,theta]=meshgrid(r,theta);
% x=r.*cos(theta);
% y=r.*sin(theta);
% z=r./tan(coneang);
x=[0,R*tan(coneang),-R*tan(coneang),0];
y=[0,R,R,0];
for j=1:1:length(x)
   XX=Rot*[x(j);y(j)]+ p;    
   x(j)=XX(1);
   y(j)=XX(2);
end
plot(x,y)
fill(x,y,'r')

% surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% alpha(0.3)




end
% keyboard
tt=0:0.0001:2*pi;
x=zeros(length(tt),1);
y=zeros(length(tt),1);
for i=1:1:length(tt)
    x(i)=Re*cos(tt(i));
    y(i)=Re*sin(tt(i));
end
plot(x,y)
fill(x,y,'b')
alpha 0.5
axis square

plot_prop_paper
% for i=1:1:Radmodel.Nsat
%     mu=XsigSat{i,1}(k,1:2);
% % Create textarrow
% [xaff,yaff] = ds2nfu([mu(1)-0.01*mu(1),mu(1)],[mu(2)-0.01*mu(2),mu(2)]);
% if max(abs(xaff))>1 || max(abs(yaff))>1
%     keyboard
% end
% annotation('textbox',[xaff(1),yaff(1),0.001,0.001],'String',{strcat('Sat ',num2str(i))});
% end

hold off


end




% 
% function txt =datatipupdatefunc(empt,event_obj,i)
% 
% pos = get(event_obj,'Position');
% txt = {['Sat ',num2str(i)]};
% 
%    
% end







