function exporttoworkspace(~,~,app)
% Export to workspace callback
prompt = {'Enter variable name:'};
title = 'Export to workspace';
lines = 1;
def = {'DEMc'};
answer = inputdlg(prompt, title, lines, def);
if isvarname(answer{1})
    assignin('base',answer{1},app.DEM);
end
end