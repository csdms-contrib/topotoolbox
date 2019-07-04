function [classname,itemname,ix] = getobjectitem(app)
% GETOBJECTITEM is used within TOPOAPP to access information about topoapp 
% objects
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013
h = findobj('Tag','Objectlist');
ix = get(h,'Value');
classname = get(h,'String');
if ~iscell(classname); classname = {classname}; end
classname = classname{ix};
h = findobj('Tag','Itemlist');
ix = get(h,'Value');
itemname = get(h,'String');
if isempty(itemname)
    itemname = '';
    ix = [];
else
    itemname = itemname(ix);
    [~,ix,~] = intersect(app.objects.(classname).names,itemname,'stable');
end
end %