function hGUI = preparegui(H1)

% prepare GUI for TTLEM
%
% Syntax
%
%     hGui = preparegui(H1)
%
%


% create figure with panels
hGUI.fig     = figure('Name','TTLEM',...
                    'NumberTitle','off',...
                    'Toolbar','none',...
                    'MenuBar','none',...
                    'Visible','off');

hGUI.main = uipanel('BackgroundColor','white',...
                'Position',[.25 .25 .75 .75],...
                'BorderType','none');

hGUI.ax.main = axes('Units','normalized', ...
               'Position', [0 0 1 1], ...
               'Parent', hGUI.main, ...
               'Clipping', 'off', ...
               'ActivePositionProperty','position',...
               'XLimMode','auto',...
               'YLimMode','auto',... 
               'PlotBoxAspectRatio',[1 1 1],...
               'PlotBoxAspectRatioMode','manual',...
               'DataAspectRatio',[1 1 1],...
               'DataAspectRatioMode','manual');%,...
daspect(hGUI.ax.main,[1 1 1]);   
hGUI.hIm = imagesc(H1.Z);  

% hold on
% xl = xlim(hGUI.ax.main);
% yl = ylim(hGUI.ax.main);


hGUI.main.SizeChangedFcn = @(varargin) FigSizeHasChanged; %set(hGUI.ax.main,'Position',[0 0 1 1]);
axis(hGUI.ax.main,'off')


% control panel
hGUI.control    = uipanel('Position',[0 0.25 0.25 0.75],...
                          'BorderType','none');
hGUI.timedisp   = uicontrol('Style', 'text', 'String', '0',...
    'units','normalized',...
    'Position', [0 .8 1 .1],...
    'Parent',hGUI.control);
hGUI.stopbutton = uicontrol('Style', 'togglebutton', 'String', 'Stop',...
    'units','normalized',...
    'Position', [0 .7 1 .1],...
    'Min',false,'Max',true,'value',false,...
    'Parent',hGUI.control);




hGUI.time     = uipanel('Position',[0 0 1 0.25],...
                        'BackgroundColor','white',...
                        'BorderType','line');
hGUI.ax.time  = axes('Parent',hGUI.time);
hold(hGUI.ax.time,'on');
xy            = get(hGUI.ax.main,'clim');
hGUI.tsmin    = plot(hGUI.ax.time,0,xy(1),'o-','Color',[.5 .5 .5],'MarkerSize',3);
hGUI.tsmean   = plot(hGUI.ax.time,0,mean(H1.Z(~isnan(H1.Z(:)))),'ko-','MarkerSize',3);
hGUI.tsmax    = plot(hGUI.ax.time,0,xy(2),'o-','Color',[.5 .5 .5],'MarkerSize',3);
legend([hGUI.tsmin; hGUI.tsmean],'min/max elevation','mean elevation');

xlabel('Time [y]');
ylabel('Elevation [m]');
set(hGUI.ax.time,'FontSize',8);

box on;

hGUI.fig.Visible = 'on';


function FigSizeHasChanged
daspect(hGUI.ax.main,[1 1 1]);
% xl = xlim(hGUI.ax.main);
% yl = ylim(hGUI.ax.main);



    
end
end
