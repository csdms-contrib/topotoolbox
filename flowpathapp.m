function flowpathapp(FD,DEM,S)

%FLOWPATHAPP Map, visualize and export flowpaths that start at manually set channelheads
%
% Syntax
%     
%     flowpathapp(FD,DEM)
%     flowpathapp(FD,DEM,S)
%
% Description
%
%     flowpathapp provides an interactive tool to visualize and generate
%     flow paths on a digital elevation model based on single flow direction
%     (FLOWobj). If an instance of STREAMobj derived from FD is supplied as
%     third argument, channelheads are automatically snapped to the
%     existing stream network, so that a subset of the latter can be
%     generated.
%
%     The stream network constructed by manually setting channelheads can
%     be exported to the workspace as new instance of STREAMobj. In
%     addition, the stream network can be exported to Excel, as text file
%     or as shapefile (requires Mapping Toolbox).
%
%     Tools are found in the menu bar of the main window. 
%
%
% Notes on export to .xls or .txt
%
%     The mapped stream network can be exported as xls or txt file. Note
%     that individual streams are listed in the order at which they were
%     mapped. Streams that are tributary to a previously mapped stream
%     terminate at the confluence of both so it may make sense to map
%     higher order stream first before mapping their tributaries. The data
%     is exported in a columnar format that include following fields:
%
%                    ID    unique stream identifier.
%                     X    x coordinates of the stream vertices
%     	              Y    y coordinates
%     upstream_distance   distance from the outlet
%             elevation   elevation of the stream vertices
%          upslope_area   drainage area of each stream vertice in map units 
%                         (nr of upstream cells x cellsize^2)
%                 slope   channel gradient in tangens [m/m]. Note that the
%                         slope is calculated along the drainage network
%                         and may differ from the slope returned by the
%                         function gradient8. If the DEM was not
%                         hydrologically conditioned, negative slope values
%                         (upward directed) may occur. 
%
%
% Input arguments
%
%     FD     FLOWobj
%     DEM    Digital elevation model (GRIDobj)
%     S      STREAMobj derived from FD
%
% 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013

% Default line color
props.linecolor = 'k';


% If you set fastindexing to true then downstream processing is much faster
% for tools that do computations along single flow paths such as
% flowpathextract
FD.fastindexing = true;

% create figure
scrsz = get(0,'ScreenSize');
% pos = [left, bottom, width, height]
hFig = figure('OuterPosition',[1/4*scrsz(3) 1/3*scrsz(4) 3/4*scrsz(3) 2/3*scrsz(4)],...
              'MenuBar','none',...
              'NumberTitle','off',...
              'Name','Main');
          
% create menu
hMenuView    = uimenu(hFig,'Label','View'); 
hButtonMenuColor = uimenu(hMenuView,'Label','Line Color'); 
colors = 'kbgyr';
hButtonColors(1) = uimenu(hButtonMenuColor,'Label','black','Checked','on','Callback',@(src,event) changelinecolor(src,event)); 
hButtonColors(2) = uimenu(hButtonMenuColor,'Label','blue','Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonColors(3) = uimenu(hButtonMenuColor,'Label','green','Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonColors(4) = uimenu(hButtonMenuColor,'Label','yellow','Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonColors(5) = uimenu(hButtonMenuColor,'Label','red','Checked','off','Callback',@(src,event) changelinecolor(src,event));


hButtonClear = uimenu(hMenuView,'Label','Clear','Callback',@(src,event) clearvectorplots);  

hMenuExport   = uimenu(hFig,'Label','Export'); 
hButtonExport = uimenu(hMenuExport,'Label','Export STREAMobj to workspace','Callback',@(src,event) exporttoworkspace);    
hButtonExportXLS = uimenu(hMenuExport,'Label','Export streams to Excel','Callback',@(src,event) writetoexcel); 
hButtonExportTXT = uimenu(hMenuExport,'Label','Export streams to ASCII','Callback',@(src,event) writetotxtfile);
hButtonExportSHP = uimenu(hMenuExport,'Label','Export streams to Shapefile','Callback',@(src,event) writetoshape);

% calculate hillshade

hs = true;
if hs;
    RGB  = imageschs(DEM,DEM,'colormap',[1 1 1]);
else
    RGB  = DEM.Z;
end
% create empty axes in figure
% hAx = axes('parent',hFig);
hAx = imgca(hFig);

if nargin == 3;
    WW = false(DEM.size);
    WW(S.IXgrid) = true;
    [rr,cc] = ind2sub(DEM.size,S.IXgrid);
    [~,~,rr,cc] = STREAMobj2XY(S,rr,cc);
%     RGB(repmat(WW,[1 1 3])) = 150;    
    [~,SNAPRASTER] = bwdist(WW,'q');
    
    FD.ixcix(~WW) = 0;
    snap = true;
else
    snap = false;
end

% show hillshade using imshow
hIm = image(RGB,'parent',hAx);
% create an instance of imscrollpanel and use the API
hPanel = imscrollpanel(hFig,hIm);
api    = iptgetapi(hPanel);
% set to initial magnification
mag    = api.findFitMag();
api.setMagnification(mag);

% create magnification box in another window
hMagBox = immagbox(hFig,hIm);
pos     = get(hMagBox,'Position');
set(hMagBox,'Position',[0 0 pos(3) pos(4)])
imoverview(hIm)

% create figure for profiles
[hFigProfiles,hAxProfiles] = getprofilefig;

if nargin == 3
    hold(hAx,'on')
    plot(hAx,cc,rr,'color',[.7 .7 .7]);
    hold(hAx,'off')
end


%% add callbacks
enterFcn = @(figHandle, currentPoint)...
       set(hFig, 'Pointer','crosshair');
iptSetPointerBehavior(hIm,enterFcn);
iptPointerManager(hFig);

xlim = get(hIm,'XData');
ylim = get(hIm,'YData');
X = 1:xlim(2);
Y = 1:ylim(2);

set(hFig,'WindowButtonDownFcn',@PressLeftButton);

%% Predefine variables
IXchannelhead = 0;
counter = 0;
LOGgrid = zeros(DEM.size,'uint16');
IXchannel = {};
distance  = {};
hPlot = [];
hPlotProfiles = [];



    function PressLeftButton(src,evt)
        
        % get the axis position
        p = get(hAx,'CurrentPoint');
        
        % check if current point is located inside axis
        p = p(1,1:2);
        
        pixelx = round(axes2pix(xlim(2), xlim, p(1)));
        pixely = round(axes2pix(ylim(2), ylim, p(2)));
        
        if pixelx < xlim(1) || pixelx > xlim(2) ...
                || pixely < ylim(1) || pixely > ylim(2)
            % do nothing
        else
            % call function
            IXchannelhead(:) = sub2ind(DEM.size,pixely,pixelx);
            counter = counter + 1;
            
            % snap to stream
            if snap
                IXchannelhead(:) = SNAPRASTER(IXchannelhead);
            end
            
            PlotProfiles
        end
        
        
    end



    function PlotProfiles
        % get coordinates of mouse click
        
        % get indices of flow path
        [IXchannel{counter},distance{counter}] = flowpathextract(FD,IXchannelhead);
        % set values in the ixcix raster (see FD.fastindexing) to zero
        % where the channel has been identified. This will ensure that
        % these locations are not visited again and thus for fast execution
        if LOGgrid(IXchannel{counter}(end)) == 0;
            LOGgrid(IXchannel{counter}) = counter;
            distance{counter} = distance{counter}(end) - distance{counter};
        else
            tribIX = LOGgrid(IXchannel{counter}(end));
            LOGgrid(IXchannel{counter}(1:end-1)) = counter;
            I = IXchannel{tribIX} == IXchannel{counter}(end);
            distance{counter} = distance{tribIX}(I)+(distance{counter}(end)-distance{counter});
        end
        
        FD.ixcix(IXchannel{counter}) = 0;
        % plot flow path
        hold(hAx,'on')
        [r,c] = ind2sub(DEM.size,IXchannel{counter});
        hPlot(counter) = plot(hAx,X(c),Y(r),props.linecolor,'LineWidth',2);
        hold(hAx,'off');
        drawnow

        hold(hAxProfiles,'on')

        hPlotProfiles(counter) = plot(hAxProfiles,distance{counter},DEM.Z(IXchannel{counter}),props.linecolor);
        hold(hAxProfiles,'off');
        drawnow
        
    end

    function clearvectorplots
        delete(hPlot);
        hPlot = [];
        delete(hPlotProfiles)
        hPlotProfiles = [];
        counter = 0;
        FD.fastindexing = true;
        
        if snap
            FD.ixcix(~WW) = 0;
        end
        
        LOGgrid = zeros(DEM.size,'uint16');
        IXchannel = {};
        distance  = {};
        
    end

    function exporttoworkspace
        W = DEM;
        W.Z = LOGgrid>0;
        
        if any(W.Z(:))
            S = STREAMobj(FD,W);
            
            prompt = {'Enter variable name:'};
            title = 'Export';
            lines = 1;
            def = {'S'};
            answer = inputdlg(prompt, title, lines, def);
            if ~isempty(answer) && isvarname(answer{1})
                assignin('base',answer{1},S);
            else
                return
            end
        else
            warndlg('No streams available for export.');
        end

    end

    function writetoexcel(src,event)
        
        
        if isempty(IXchannel)
            warndlg('No streams available for export.');
        else
            [D,header] = makedataset(IXchannel,distance);
            D = [header; num2cell(D)];
            
            [FileName,PathName] = uiputfile({'*.xlsx';'*.xls'},'Write to Excel');
            
            if FileName == 0
                return
            end
            
            xlswrite([PathName FileName],D);
        
        end
    end

    function writetotxtfile(src,event)
 
        if isempty(IXchannel)
            warndlg('No streams available for export.');
           
        else
            [D,header] = makedataset(IXchannel,distance);            
            [FileName,PathName] = uiputfile({'*.txt'},'Write to text file');
            
            if FileName == 0
                return
            end
            
            fid = fopen([PathName FileName], 'w');
            
            for r = 1:numel(header);
                fprintf(fid, header{r});
                if r < numel(header)
                    fprintf(fid, '\t');
                end 
            end
            fprintf(fid, '\n');

            for row=1:size(D,1);
                fprintf(fid, '%d\t%f\t%f\t%f\t%f\t%f\t%f\n', D(row,:));
            end

            fclose(fid);
        
        end
    end

    function writetoshape(src,event)
        if isempty(IXchannel)
            warndlg('No streams available for export.');
        elseif ~(exist('shapewrite','file')==2);
            warndlg('The function shapewrite is not available. Shapewrite is part of the Mapping Toolbox.');
        else
            
            
            for r = 1:numel(IXchannel);
                SHP(r).Geometry = 'Line';
                [SHP(r).X SHP(r).Y] = ind2coord(DEM,IXchannel{r}(:)); 
                SHP(r).X = SHP(r).X';
                SHP(r).Y = SHP(r).Y';
                SHP(r).ID = r;
                SHP(r).minZ = double(min(DEM.Z(IXchannel{r})));
                SHP(r).maxZ = double(max(DEM.Z(IXchannel{r})));
                SHP(r).length = double(max(distance{r})-min(distance{r}));
                SHP(r).tribtoID = double(LOGgrid(IXchannel{r}(end)));
                if SHP(r).tribtoID == SHP(r).ID;
                    SHP(r).tribtoID = 0;
                end
            end
            [FileName,PathName] = uiputfile('*.shp','Write to Shapefile');
            
            if FileName == 0
                return
            end
            
            shapewrite(SHP,[PathName FileName]);
        end
  
    end

    function changelinecolor(src,event)
        
        
        I = src == hButtonColors;

        props.linecolor = colors(I);

        
        h = findobj(hButtonColors,'Checked','on');
        set(h,'Checked','off');
        set(src,'Checked','on');
    end
  

    function [D,header] = makedataset(IXchannel,distance)
                
        A = flowacc(FD);
        G = FLOWobj2gradient(FD,DEM);
        D = cell(numel(IXchannel),1);
        for r = 1:numel(IXchannel);
            D{r}(1:numel(IXchannel{r}),1) = repmat(r,numel(IXchannel{r}),1);
            [D{r}(:,2) D{r}(:,3)] = ind2coord(DEM,IXchannel{r}(:));            
            D{r}(:,4) = distance{r}(:);
            D{r}(:,5) = DEM.Z(IXchannel{r}(:));
            D{r}(:,6) = A.Z(IXchannel{r}(:)).*(DEM.cellsize).^2;
            D{r}(:,7) = G.Z(IXchannel{r}(:));
        end
        
        D = cell2mat(D);
        header = {'ID' 'X' 'Y' 'upstream_distance' 'elevation' 'upslope_area' 'slope'};
   
    end
    
    function [hFigProfiles,hAxProfiles] = getprofilefig
        % create figure for profiles
        hFigProfiles = figure('OuterPosition',[1/4*scrsz(3) 50 3/4*scrsz(3) 1/3*scrsz(4)-50],...
                      'NumberTitle','off',...
                      'Name','Profiles') ;
        hAxProfiles = axes('Parent',hFigProfiles,'Xscale','linear','Yscale','linear','box','on');
        xlabel(hAxProfiles,'distance from outlet');
        ylabel(hAxProfiles,'elevation');
    end    

        
end
        





