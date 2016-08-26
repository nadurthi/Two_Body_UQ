function doc_datacursormode
% Plots graph and sets up a custom data tip update function
fig = figure;
a = -16; t = 0:60;
plot(t,sin(a*t))
dcm_obj = datacursormode(fig);
h = createDatatip(dcm_obj,fig);
drawnow
set(h,'Position',[1 2],...
  'String',sprintf('X: %s\nY: %s',...
  num2str(1),num2str(2)));

% set(dcm_obj,'UpdateFcn',@myupdatefcn)
% datacursormode on

end


