function enabletools(app,state)

% toggle availability of tools
%
% enabletools(app,state)
%
% State can be  'off' or 'on'
set(app.menu.methods,'enable',state)
set(app.menu.view,'enable',state)
set(app.menu.export,'enable',state)
set(app.menu.edit,'enable',state)
set(app.menu.calc,'enable',state)
set(app.tbhb(:),'enable',state)
zoom(app.fh,state)
pan(app.fh,state)
pause(.1)
end