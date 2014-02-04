function layout(app)

% create layout for gui_preprocess
%
% layout(app)
%
% 

%% set up figure
app.fh = figure('units','pixels',...
    'position',[100 100 800 600],...
    'menubar','none',...
    'name','DEM preprocessing tool',...
    'numbertitle','off',...
    'resize','on',...
    'ResizeFcn',{@figresize, app},...
    'toolbar','figure',...
    'CloseRequestFcn',{@closefcn, app});

%% set up menu
app.menu.methods    = uimenu(app.fh,'Label','Methods');
app.menu.menucarve  = uimenu(app.menu.methods,'Label','Carving','Checked','on','Callback',{@showtools, app});
app.menu.menufill   = uimenu(app.menu.methods,'Label','Filling','Checked','off','Callback',{@showtools, app});

app.menu.edit       = uimenu(app.fh,'Label','Edit');
app.menu.menuundo   = uimenu(app.menu.edit,'Label','Undo','Enable','off','Callback',{@undo, app});

app.menu.view       = uimenu(app.fh,'Label','View');
app.menu.menuss     = uimenu(app.menu.view,'Label','Show sinks ','Checked','off','Callback',{@showsinks, app});
app.menu.menuim     = uimenu(app.menu.view,'Label','Show scaled image ','Checked','on','Separator','on','Callback',{@showimage, app});
app.menu.menuac     = uimenu(app.menu.view,'Label','Adjust color range','Checked','off','Callback',{@adjustcolorrange, app});
app.menu.menuhs     = uimenu(app.menu.view,'Label','Show hillshade ','Checked','off','Separator','on','Callback',{@showhillshade, app});

app.menu.calc       = uimenu(app.fh,'Label','Calculate');
app.menu.menucalcsinks  = uimenu(app.menu.calc,'Label','Recalculate sinks','Callback',{@calculatesinks, app});
app.menu.menucalchs     = uimenu(app.menu.calc,'Label','Recalculate hillshade','Callback',{@calculatehs, app});

app.menu.export     = uimenu(app.fh,'Label','Export');
app.menu.menuexport = uimenu(app.menu.export,'Label','Export DEM to workspace','Callback',{@exporttoworkspace, app});


%%
% set up zoom toolbar
tbh      = findall(app.fh,'Type','uitoolbar');
app.tbhb = findall(tbh);
delete(app.tbhb([2:9 13:end]));
app.tbhb = app.tbhb(10:12);

% axes
app.ax = axes('parent',app.fh,...
    'Units','pixels',...
    'DataAspectRatio',[1 1 1],...
    'DataAspectRatioMode','manual',...
    'Position',[0 0 600 600],...
    'Xlimmode','auto',...
    'Ylimmode','auto',...
    'PlotBoxAspectRatio',[1 1 1],...
    'PlotBoxAspectRatioMode','manual',...
    'Visible','off',...
    'Layer','top',...
    'Box','on');


% display DEM in axes
app.im = image('parent',app.ax,...
    'XData',app.X,...
    'YData',app.Y,...
    'CData',app.DEM.Z,...
    'CDataMapping','scaled',...
    'visible','on');

% zoom options
%Setup zoom
set(app.ax,'ylimmode','auto','xlimmode','auto'); %reset after image detroyed
app.zm = zoom(app.ax); %build zoom object
zoom reset; %store current setting
set(app.zm,'ActionPostCallback',{@mypostcallback,app}); %set callback

% set up buttons
app.mpanel = uipanel('Units','Pixels',...
    'Position',[601 1 199 598]);

% text field for help text
app.helptext = uicontrol('Parent',app.mpanel,...
    'Style','text',...
    'Unit','normalized',...
    'Position',[0.02 0.91 0.96 0.08],...
    'String','Choose between carving and filling',...
    'HorizontalAlignment','left',...
    'BackgroundColor',[.9 .9 .9]);

%% ***********************************************************************
% toggle buttons for switching between carving and filling
app.methodgroup = uibuttongroup('Parent',app.mpanel,...
    'Position',[0 0.85 1 0.05]);

%% ***********************************************************************
% panel for carving
app.panelcarve = uipanel('Parent',app.mpanel,...
    'Unit','normalized',...
    'Position',[0 0.6 1 0.3],...
    'visible','on');
app.carvegroup = uibuttongroup('Parent',app.panelcarve,...
    'Unit','normalized',...
    'Position',[0 0.5 1 0.5],...
    'bordertype','none');

app.carvemeth(1) = uicontrol('Style','Radio','String','linear interpolation','Unit','normalized',...
    'pos',[0 0.66+0.1 1 0.2],'parent',app.carvegroup,'HandleVisibility','off');
app.carvemeth(2) = uicontrol('Style','Radio','String','linear minima imposition','Unit','normalized',...
    'pos',[0 0.33+0.1 1 0.2],'parent',app.carvegroup,'HandleVisibility','off');
app.carvemeth(3) = uicontrol('Style','Radio','String','least cost path minima imposition','Unit','normalized',...
    'pos',[0 0+0.1 1 0.2],'parent',app.carvegroup,'HandleVisibility','off');

app.carvevaltext = uicontrol('Parent',app.panelcarve,...
    'Style','text',...
    'Unit','normalized',...
    'Position',[0 0.2 0.5 0.2],...
    'String','rel. box size = ',...
    'HorizontalAlignment','right');
app.carvevaledit = uicontrol('Parent',app.panelcarve,...
    'Style','edit',...
    'Unit','normalized',...
    'Position',[0.55 0.28 0.3 0.13],...
    'String','1.2',...
    'HorizontalAlignment','right',...
    'backgroundcolor','white');

% toggle button for carving
app.bcarve = uicontrol('parent',app.panelcarve',...
    'Style','togglebutton',...
    'Unit','normalized',...
    'position',[0 0 1 0.2],...
    'string','Start',...
    'Callback',{@startcarving, app});

%% ***********************************************************************
% panel for filling
app.panelfill = uipanel('Parent',app.mpanel,...
    'Unit','normalized',...
    'Position',[0 0.6 1 0.3],...
    'visible','off');
app.fillgroup = uibuttongroup('Parent',app.panelfill,...
    'Unit','normalized',...
    'Position',[0 0.5 1 0.5],...
    'bordertype','none');

app.fillmeth(1) = uicontrol('Style','Radio','String','Fill sinks (areas less than [px])','Unit','normalized',...
    'pos',[0 0.66+0.1 1 0.2],'parent',app.fillgroup,'HandleVisibility','off');
app.fillmeth(2) = uicontrol('Style','Radio','String','Fill sinks (depth less than [z units])','Unit','normalized',...
    'pos',[0 0.33+0.1 1 0.2],'parent',app.fillgroup,'HandleVisibility','off');
app.fillmeth(3) = uicontrol('Style','Radio','String','Fill sinks manually','Unit','normalized',...
    'pos',[0 0+0.1 1 0.2],'parent',app.fillgroup,'HandleVisibility','off');

app.fillvaltext = uicontrol('Parent',app.panelfill,...
    'Style','text',...
    'Unit','normalized',...
    'Position',[0 0.2 0.5 0.2],...
    'String','area or depth = ',...
    'HorizontalAlignment','right');
app.fillvaledit = uicontrol('Parent',app.panelfill,...
    'Style','edit',...
    'Unit','normalized',...
    'Position',[0.55 0.28 0.3 0.13],...
    'String','1',...
    'HorizontalAlignment','right',...
    'backgroundcolor','white');

% toggle button for filling
app.bfill = uicontrol('parent',app.panelfill',...
    'Style','togglebutton',...
    'Unit','normalized',...
    'position',[0 0 1 0.2],...
    'string','Start',...
    'Callback',{@startfilling, app});


end