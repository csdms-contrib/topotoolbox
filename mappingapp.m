function mappingapp(DEM,S,varargin)

% map knickpoints combining planform and profile view
%
% Syntax
%
%     mappingapp(DEM,S)
%
% Description
%
%     This light-weight tool enables mapping points in planform and profile
%     view simultaneously based on a digital elevation model (DEM) and a
%     stream network (S). Points are stored as a table and can be exported.
%     Note that this tool is still beta and quite rudimentary.
%
% Input arguments
%
%     DEM    Digital elevation model (GRIDobj)
%     S      Stream network (STREAMobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     mappingapp(DEM,S)    
%
%
% See also: STREAMobj, flowpathapop
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 7. November, 2018


% Add the GUI components
hs = addcomponents;

% Add data to GUI
% calculate hillshade
RGB  = imageschs(DEM,DEM,'colormap',landcolor(255));
% create empty axes in figure
ax.main = axes('Units','normalized', ...
               'Position', [0 0 1 1], ...
               'Parent', hs.main, ...
               'Clipping', 'off', ...
               'ActivePositionProperty','position');

% After adding all components make figure visible
hs.fig.Visible = 'on';           
% % set axis to image
% axis(ax.main,'image');
% show hillshade using imshow
hIm = imshow(RGB,'parent',ax.main);

% create an instance of imscrollpanel and use the API
hPanel = imscrollpanel(hs.main,hIm);
set(hPanel,'Units','normalized','Position',[0 0 1 1])
api    = iptgetapi(hPanel);

hs.Menu.ZoomIn.ClickedCallback = @callbackZoomIn;
hs.Menu.ZoomOut.ClickedCallback = @callbackZoomOut;

hImOverview = imoverviewpanel(hs.overview,hIm);           
           
% get s
[x,y,d,z] = STREAMobj2XY(S,S.distance,DEM);
[r,c] = coord2sub(DEM,x,y);

hold(ax.main,'on')
ls.main = plot(ax.main,c,r,'k');
hold(ax.main,'off')

% profile
ax.profile = axes('parent',hs.profiles,...
    'Position',[0 0 1 1],...
    'TickLength',[0.005 0.005]);
ls.profile = plot(ax.profile,d,z,'Color',[.2 .2 .2]);

ix      = [];
point   = [];
activerow   = 1;
activecolor   = 'b';
inactivecolor = 'w';

% set the first row
hs.table.Data = cell(1,4);
setnewpoints;

hs.Menu.Add.ClickedCallback    = @callbackAdd;
hs.Menu.Export.ClickedCallback = @callbackExport;

%% Table utilities
function setnewpoints(pos)
    if nargin == 0
        point.main = impoint(ax.main,'PositionConstraintFcn',@getnearestmain);
    else
        point.main = impoint(ax.main,pos,'PositionConstraintFcn',@getnearestmain);
    end
    setColor(point.main,activecolor)
    setPosition(point.main,getPosition(point.main))

    point.profile = impoint(ax.profile,[d(ix) z(ix)],'PositionConstraintFcn',@getnearestprofile);
    setColor(point.profile,activecolor)
end

function callbackAdd(varargin)    
    % Add row
    hs.table.Data = [hs.table.Data;cell(1,4)];
    activerow = activerow + 1;
    hold(ax.profile,'on');
    data = hs.table.Data;
    nrEntries = size(data,1);
    inactiverows = setdiff(1:nrEntries,activerow);
    [~,I] = ismember(cell2mat(data(inactiverows,1:2)),[x y],'rows');
    hold(ax.main,'on');
    plot(ax.main,c(I),r(I),'ok','MarkerFaceColor',inactivecolor);
    hold(ax.main,'off');
    
    hold(ax.profile,'on');
    plot(ax.profile,d(I),z(I),'ok','MarkerFaceColor',inactivecolor);
    hold(ax.profile,'off');
    
    % Now move the mobile point a bit downstream because otherwise it is
    % difficult to grab
    setnewpoints(getPosition(point.main))
%     setPosition(point.main,getPosition(point.main))
%     setPosition(point.profile,getPosition(point.profile));
end



function callbackZoomIn(varargin)
    % Zoom in
    api.setMagnification(api.getMagnification() + 1);
end

function callbackZoomOut(varargin)
    % Zoom out
    mag = api.findFitMag();
    api.setMagnification(max(api.getMagnification() - 1,mag));
end

function callbackExport(varargin)
    % Create table
    T = cell2table(hs.table.Data,'VariableNames',hs.table.ColumnName);
    
    if iscell(T.X(end))
        T = cell2table(hs.table.Data(1:end-1,:),'VariableNames',hs.table.ColumnName);
    end
    
    % get linear indices
    [~,locb] = ismember([T.X T.Y],[S.x S.y],'rows');
    T.IXgrid = S.IXgrid(locb); 

    prompt = {'Enter variable name:'};
    title = 'Export';
    lines = 1;
    def = {'T'};
    answer = inputdlg(prompt, title, lines, def);
    if ~isempty(answer) && isvarname(answer{1})
        assignin('base',answer{1},T);
    else
        return
    end

end

function posn = getnearestmain(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((c-pos(1)).^2 + (r-pos(2)).^2);
    posn= [c(ix) r(ix)];    
    if isfield(point,'profile') 
        setPosition(point.profile,[d(ix) z(ix)]);
    end    
    hs.table.Data{activerow,1} = x(ix);
    hs.table.Data{activerow,2} = y(ix);
    hs.table.Data{activerow,3} = z(ix);
end

function posn = getnearestprofile(pos)
    % GETNEAREST  snap to nearest location on stream network profile
    [~,ix] = min((d-pos(1)).^2 + (z-pos(2)).^2);
    posn= [d(ix) z(ix)];
    setPosition(point.main,[c(ix) r(ix)]);
    hs.table.Data{activerow,1} = x(ix);
    hs.table.Data{activerow,2} = y(ix);
    hs.table.Data{activerow,3} = z(ix);
end

end


function hs = addcomponents
% create figure with panels
hs.fig     = figure('Name','MappingApp',...
                    'NumberTitle','off',...
                    'Toolbar','none',...
                    'MenuBar','none',...
                    'Visible','off');

hs.toolbar = uitoolbar(hs.fig);
iconin = createIcon('in');
hs.Menu.ZoomIn = uipushtool(hs.toolbar,...
    'TooltipString','Zoom In',...
    'Cdata',iconin);

iconout = createIcon('out');
hs.Menu.ZoomOut = uipushtool(hs.toolbar,...
    'TooltipString','Zoom Out',...
    'Cdata',iconout);

iconadd = createIcon('add');
hs.Menu.Add = uipushtool(hs.toolbar,...
    'TooltipString','Add row',...
    'Cdata',iconadd,...
    'Separator','on');

% icondel = createIcon('del');
% hs.Menu.Del = uipushtool(hs.toolbar,...
%     'TooltipString','Delete row',...
%     'Cdata',icondel);


iconexport = createIcon('export');
hs.Menu.Export = uipushtool(hs.toolbar,...
    'TooltipString','Export to workspace',...
    'Cdata',iconexport);


% main panel with hillshade
hs.main = uipanel('BackgroundColor','white',...
                'Position',[.25 .25 .75 .75]);
            
% overview panel to help navigation
hs.overview = uipanel('Position',[0 0.75 0.25 0.25]); 

% table panel
hs.features = uipanel('Position',[0 0.25 0.25 0.50]);

% Column names and column format
fcolumnname = {'X','Y','Z','Name'};
fcolumnformat = {'numeric','numeric','numeric','char'};
hs.table = uitable('Parent',hs.features,...
            'ColumnName', fcolumnname,...
            'ColumnFormat', fcolumnformat,...
            'ColumnEditable', [false false false true],...
            'Units','normalize',...
            'Position',[0 0 1 1]);
        
% profile panel
hs.profiles = uipanel('Position',[0 0 1 .25]);

end
         

function c = createIcon(type)

switch type
    case 'in'  
        c = [...
                0    0    2   84  173  215  214  170   78    1    0    0    0    0    0    0
                0   18  187  252  253  253  253  253  252  178   14    0    0    0    0    0
                5  193  253  248  154   87   89  160  250  253  182    2    0    0    0    0
              102  252  245   65    3  181  181    3   75  248  252   88    0    0    0    0
              198  253  135    0    3  237  237    3    0  148  253  182    1    0    0    0
              245  253   58  179  210  251  251  210  179   71  253  229    2    0    0    0
              249  253   51  228  252  253  253  252  229   63  253  233    2    0    0    0
              212  253  111    8   16  237  237   16    8  125  253  196    1    0    0    0
              127  253  230   30    0  219  219    0   37  236  253  114    2    0    0    0
               17  222  253  226   98   45   47  104  231  254  254  209  151   14    0    0
                0   45  225  253  253  253  253  253  253  254  254  254  253  194   14    0
                0    0   20  139  227  252  252  224  134  183  254  255  254  253  194   15
                0    0    0    1    2   18   17    1    2  100  253  254  255  254  253  175
                0    0    0    0    0    0    0    0    0    1  139  253  254  255  254  249
                0    0    0    0    0    0    0    0    0    0    1  139  253  254  253  196
                0    0    0    0    0    0    0    0    0    0    0    1  118  213  177   31];
    case 'out'
        c = [...
                0    0    2   84  173  215  214  170   78    1    0    0    0    0    0    0
                0   18  187  252  253  253  253  253  252  178   14    0    0    0    0    0
                5  193  253  248  153   86   87  159  250  253  182    2    0    0    0    0
              102  252  245   65    1    0    0    1   75  248  252   88    0    0    0    0
              198  253  135    2    2    1    1    2    2  148  253  182    1    0    0    0
              245  253   58  179  209  209  209  209  179   71  253  229    2    0    0    0
              249  253   51  228  252  252  252  252  229   63  253  233    2    0    0    0
              212  253  111    8   15   15   15   15    8  125  253  196    1    0    0    0
              127  253  230   30    1    0    0    1   37  236  253  114    2    0    0    0
               17  222  253  226   97   30   32  103  231  254  254  209  151   14    0    0
                0   45  225  253  253  253  253  253  253  254  254  254  253  194   14    0
                0    0   20  139  227  252  252  224  134  183  254  255  254  253  194   15
                0    0    0    1    2   18   17    1    2  100  253  254  255  254  253  175
                0    0    0    0    0    0    0    0    0    1  139  253  254  255  254  249
                0    0    0    0    0    0    0    0    0    0    1  139  253  254  253  196
                0    0    0    0    0    0    0    0    0    0    0    1  118  213  177   31];
    case 'add'
        c = [ ...
                0    0    0    0    0   12  178  251  251  178   12    0    0    0    0    0
                0    0    0    0    0  105  253  254  254  253  105    0    0    0    0    0
                0    0    0    0    1  126  253  255  255  253  126    1    0    0    0    0
                0    0    0    0    1  127  253  255  255  253  127    1    0    0    0    0
                0    0    1    1    2  127  254  255  255  254  127    2    1    1    0    0
               12  105  126  127  127  191  254  255  255  254  191  127  127  126  105   12
              178  253  253  253  254  254  254  255  255  254  254  254  253  253  253  178
              251  254  255  255  255  255  255  255  255  255  255  255  255  255  254  251
              251  254  255  255  255  255  255  255  255  255  255  255  255  255  254  251
              178  253  253  253  254  254  254  255  255  254  254  254  253  253  253  178
               12  105  126  127  127  191  254  255  255  254  191  127  127  126  105   12
                0    0    1    1    2  127  254  255  255  254  127    2    1    1    0    0
                0    0    0    0    1  127  253  255  255  253  127    1    0    0    0    0
                0    0    0    0    1  126  253  255  255  253  126    1    0    0    0    0
                0    0    0    0    0  105  253  254  254  253  105    0    0    0    0    0
                0    0    0    0    0   12  178  251  251  178   12    0    0    0    0    0];
        
    case 'del'
        c = [...
                1   65  217  248  166   13    0    0    0    0   13  166  248  217   66    1
               65  243  254  254  253  190   12    0    0   12  190  253  254  254  244   66
              217  254  254  255  254  253  190   13   13  190  253  254  255  254  254  217
              248  254  255  255  255  254  253  191  191  253  254  255  255  255  254  248
              166  253  254  255  255  255  254  254  254  254  255  255  255  254  253  166
               13  190  253  254  255  255  255  255  255  255  255  255  254  253  190   13
                0   12  190  253  254  255  255  255  255  255  255  254  253  190   12    0
                0    0   13  191  254  255  255  255  255  255  255  254  191   13    0    0
                0    0   13  191  254  255  255  255  255  255  255  254  191   13    0    0
                0   12  190  253  254  255  255  255  255  255  255  254  253  190   12    0
               13  190  253  254  255  255  255  255  255  255  255  255  254  253  190   13
              166  253  254  255  255  255  254  254  254  254  255  255  255  254  253  166
              248  254  255  255  255  254  253  191  191  253  254  255  255  255  254  248
              217  254  254  255  254  253  190   13   13  190  253  254  255  254  254  217
               66  244  254  254  253  190   12    0    0   12  190  253  254  254  244   66
                1   66  217  248  166   13    0    0    0    0   13  166  248  217   66    1];
    case 'export'
        c = [ ...
                0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
                0    0    0    0    0    0    0   96  255  255  255  255  255  255  255    0
                0    0    0    0    0    0    0    0   96  255  255  255  255  255  255    0
                0    0    0    0    0    0    0    0    0  191  255  255  255  255  255    0
               63  255  255  255  255  255  175   16  159  255  255  255  255  255  255    0
              255  255  255  255  255  239   64  159  255  255  255  255  255  255  255    0
              255  255   63    0    0   16  207  255  255  255  255  255  159  255  255    0
              255  255    0    0   16  207  255  255  255  255  255   96    0   96  255    0
              255  255    0    0   80  255  255  255  255  255   96    0    0    0   80    0
              255  255    0    0    0   96  255  255  255   96   96   64    0    0    0    0
              255  255    0    0    0    0   96  255   96   80  255  255    0    0    0    0
              255  255    0    0    0    0    0    0    0    0  255  255    0    0    0    0
              255  255    0    0    0    0    0   32    0    0  255  255    0    0    0    0
              255  255   63    0    0    0    0    0    0   63  255  255    0    0    0    0
              255  255  255  255  255  255  255  255  255  255  255  255    0    0    0    0
               63  255  255  255  255  255  255  255  255  255  255   63    0    0    0    0];

 
end
c(c<=18) = nan;
c = double(c)/255;
c(c==0) = nan;
c = 1-c;

c = repmat(c,1,1,3);
end

