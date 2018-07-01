function [h,posYlast] = listfigure(app,classes)
% LISTFIGURE creates a figure within topoapp session for access to topoapp 
% object data. Some standard fields and functionalities are included.
% Advanced functionalities can be added. posYlast gives the last position
% in the right-most column below which to add new parameters controls.
%
% See also: CHIPLOTTING, SWATHMAPPING, MANAGEOBJECTS
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

set(app.gui.TB,'Enable','off');

h = figure('Name','Topoapp: Objects','Menu','none',...
    'CloseRequestFcn',{@closelistfig,app},...
    'BusyAction','cancel');

%--------------------------------------------------------
% positioning:
dLs = 0.025; dLl = 0.2833;
Lft1 = dLs; Lft2 = Lft1+(dLl+dLs+dLs); Lft3 = Lft2+(dLl+dLs+dLs);
dLftedit1 = 0.12;  dLftedit2 = 0.2133;
Hpb = 0.06; Hedit = 0.045; Htxt = 0.045;
Wpb = dLl; Wedit = 0.07; Wtxts = 0.1; Wtxtl = 0.17;
%--------------------------------------------------------

%% UI panels
%--------------------------------------------------------
uipanel('Title','Object selection','Position',[0 0 1/3 1],...
    'Backgroundcolor',[0.8 0.8 0.8],'Fontsize',10);
uipanel('Title','Tools','Position',[1/3 0.4 1/3 0.6],...
    'Backgroundcolor',[0.8 0.8 0.8],'Fontsize',10);
uipanel('Title','Tool options','Position',[1/3 0 1/3 0.4],...
    'Backgroundcolor',[0.8 0.8 0.8],'Fontsize',10);
uipanel('Title','Parameter controls','Position',[2/3 0 1/3 1],...
    'Backgroundcolor',[0.8 0.8 0.8],'Fontsize',10);
%--------------------------------------------------------

%% Object selection
%--------------------------------------------------------
uicontrol('Style','popupmenu',...
    'units','normalized','position',[Lft1,.875,Wpb,Hpb],...
    'String',classes,...
    'Min',1,'Max',1,...
    'Tag','Objectlist',...
    'Callback',{@selectobject,app});
uicontrol('Style','listbox',...
    'units','normalized','position',[Lft1,.275,Wpb,.58],...
    'Max',999,...
    'Tag','Itemlist',...
    'Callback',{@selectitem,app},...
    'ButtonDownFcn',{@renameitem,app});
uicontrol('Style','pushbutton',...
    'units','normalized','position',[Lft1,.175,Wpb,Hpb],...
    'String','Delete object',...
    'Callback',{@deleteitem,app});
uicontrol('Style','pushbutton',...
    'units','normalized','position',[Lft1,.1,Wpb,Hpb],...
    'String','Set visible',...
    'Callback',{@showitem,app});
%--------------------------------------------------------

%% Axis controls
%--------------------------------------------------------
posY = 0.875;
uicontrol('Style','text',...
    'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
    'String','X limits','Backgroundcolor',[.8 .8 .8]);
uicontrol('Style','edit','Tag','xmin',...
    'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
    'String','min');
uicontrol('Style','edit','Tag','xmax',...
    'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
    'String','max');
posY = posY - 0.05;
uicontrol('Style','text',...
    'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
    'String','Y limits','Backgroundcolor',[.8 .8 .8]);
uicontrol('Style','edit','Tag','ymin',...
    'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
    'String','min');
uicontrol('Style','edit','Tag','ymax',...
    'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
    'String','max');

posYlast = posY - 0.05;

end

%------------------------------------------------------------
function selectobject(~,~,app)
h = findobj(gcf,'Tag','Objectlist');
classname = get(h,'String');
if ~iscell(classname); classname = {classname}; end
ix = get(h,'Value');
classname = classname{ix};
deletetempobj(app)
updatelist(classname,app)
end % 
%------------------------------------------------------------
function updatelist(classname,app)
if ~isempty(app.objects.(classname))
    ix = [app.objects.(classname).visible];
    itemnames = app.objects.(classname).names(logical(ix));
end
hitemlist = findobj(gcf,'Tag','Itemlist');
set(hitemlist,'String',itemnames,'Value',1) % display items
deselectallitems(app)
h = findobj(gcf,'Style','togglebutton');
set(h,'Value',0)
% remove active parameter controls
h = findobj(gcf,'Tag','optionalparameter');
if ~isempty(h)
    set(h,'Visible','off','Enable','off');
end
% enable optional tools
h = findobj(gcf,'Style','togglebutton','Tag','optionaltool');
if ~isempty(h)
    set(h,'Visible','off','Enable','off');
    for i = 1 : length(h)
        userdata = get(h(i),'userdata');
        if sum(strcmp(userdata.classes,classname))>0
            set(h(i),'Visible','on','Enable','on');
        end
    end
end
end % 
%------------------------------------------------------------
function closelistfig(~,~,app)
deselectallitems(app)
set(app.gui.TB,'Enable','on');
deletetempobj(app)
delete(gcf)
end % 
%------------------------------------------------------------
function selectitem(~,~,app)
hfig = gcf;
figure(app.gui.hfig)
deselectallitems(app)
[classname,~,ix] = getobjectitem(app);
h = [app.objects.(classname).handles(ix)];
set(h,'Selected','on')
figure(hfig)
end % 
%------------------------------------------------------------
function deselectallitems(app)
hc = get(app.gui.hax,'Children');
if ~isempty(hc)
    set(hc,'Selected','off')
end
end % 
%------------------------------------------------------------
function deleteitem(~,~,app)
[objectname,~,ix] = getobjectitem(app);
if ~isempty(ix)
    delete(app.objects.(objectname).handles(ix));
    app.objects.(objectname).names(ix) = [];
    app.objects.(objectname).handles(ix) = [];
    app.objects.(objectname).visible(ix) = [];
    app.objects.(objectname).data(ix) = [];
    updatelist(objectname,app)
else error('No object to delete')
end
end % 
%------------------------------------------------------------
function renameitem(~,~,app)
[objectname,itemname,ix] = getobjectitem(app);
answer = inputdlg('Enter object name:','renameitem object',1,itemname);
if ~isempty(answer)
    h = app.objects.(objectname).handles(ix);
    set(h,'DisplayName',answer{1});
    app.objects.(objectname).names{ix} = answer{1};
    updatelist(objectname,app)
end
end % 
%------------------------------------------------------------
function showitem(~,~,app)
[objectname,~,~] = getobjectitem(app);
nobjects = length(app.objects.(objectname).data);
if nobjects>0
    items = app.objects.(objectname).names;
    set(0,'DefaultUIcontrolUnits','pixel');
    [s,ok] = listdlg('PromptString','Select objects to show:',...
        'SelectionMode','multiple',...
        'InitialValue',[],...
        'ListString',items);
    set(0,'DefaultUIcontrolUnits','normalized');
    if ok
        for i = 1 : nobjects
            app.objects.(objectname).visible(i) = 0;
            set(app.objects.(objectname).handles(i),'Visible','off')
        end
        for i = 1 : length(s)
            app.objects.(objectname).visible(s(i)) = 1;
            set(app.objects.(objectname).handles(s(i)),'Visible','on')
        end
    end
    updatelist(objectname,app)
end
end % 
%------------------------------------------------------------
function deletetempobj(app)
h = findobj(app.gui.hax,'Tag','TemporaryObject');
if ~isempty(h)
    delete(h);
end
end