classdef topoapp < handle
% TOPOAPP createa an instance of class topoapp
%
% Syntax
%
%     APP = topoapp
%     APP = topoapp(DEM)
%     APP = topoapp(APP)
%
% Description
%
%     TOPOAPP creates an instance of the topoapp class, which is a 
%     graphical user interface for working with digital elevation models.
%     TOPOAPP acts as an interface for quick access to tools that are
%     offered by the various TopoToolbox functions. It also provides a
%     platform for implementing higher-order tools that comply with the
%     object classes of the TopoToolbox. Users can create their own tools
%     and append them to the functionality of the topoapp. Several tools to
%     create basic features (objects) are supplied with the TOPOAPP, such
%     as 'pickprofile', 'importoutlets', 'pickdrainage' or 'pickreach',
%     which are all located in the PRIVATE folder of the TOPOAPP. Several
%     functions, e.g., to initialize a class that is used in a tool
%     ('initclass'), to add a newly created object to the TOPOAPP
%     ('addobject'), or to create a simple figure that lists the objects
%     created in a TOPOAPP session ('listfigure'), are also supplied and
%     help to conform with the way TOPOAPP is handling the data when coding
%     new tools. The best way to learn how to implement a new tool is by
%     examining the currently available functions and tools. 
%     
%     
%     APP = topoapp      creates a topoapp object by opening a dialog for
%                        selecting a file that contains a DEM
%
%     APP = topoapp(DEM) creates a topoapp object from a DEM that must be 
%                        an instance of class GRIDobj
%
%     APP = topoapp(APP) restores a topoapp object previously created
%
%     
% 
% TOPOAPP properties:
%     
%     DEM       - digital elevation model (GRIDobj)
%     X         - coordinate vector
%     Y         - coordinate vector
%     
%     G         - slope map (GRIDobj)
%     RGB       - hillshade image colored by elevation (RGB)
%     HS        - hillshade image (grayscale)
%     HSG       - hillshade image colored by elevation slope (RGB)
%     SINKS     - hillshade image with colored sink depths (RGB)
%
%     FD        - flow direction object (FLOWobj)
%     FA        - flow accumulation grid (GRIDobj)
%     DEMc      - digital elevation model with imposed minima (GRIDobj)
%     S         - stream network (STREAMobj)
%
%     gui       - gui-related handles
%     objects   - objects created during topoapp session
%     paras     - parameter values
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


properties
    
    DEM     % digital elevation model (GRIDobj)
    X       % coordinate vector
    Y       % coordinate vector
    G       % slope map
    RGB     % hillshade colored by elevation (RGB)
    HS      % hillshade (grayscale)
    HSG     % hillshade colored by elevation slope
    SINKS   % hillshade with colored sink depths
    FD      % flow direction object (FLOWobj)
    FA      % flow accumulation grid (GRIDobj)
    DEMc    % digital elevation model with imposed minima (GRIDobj)
    S       % stream network (STREAMobj)
    gui     % structure of gui-related handles
    objects % structure of object data
    paras   % structure of parameter data
    
end % properties


methods
    
    % App methods
    %----------------------------------------------------------------
    function this_app = topoapp(varargin)
        
        % check input
        narginchk(0,1)
        if nargin==0
            [demfile,dempath] = uigetfile('.tif','Choose DEM file');
            this_app.DEM = GRIDobj([dempath,demfile]);
        elseif nargin==1
            this_app.DEM = varargin{1};
        end
        
        if isa(this_app.DEM,'topoapp')
            % restore instance of topoapp
            this_app = this_app.DEM;
            guilayout(this_app) % create graphical user interface
            
            hold on
            h = plotobject(this_app,this_app.S,'k');
            set(h,'Displayname','Base STREAMobj','visible','off');
            
            classes = this_app.objects.classes;
            for i = 1 : length(classes)
                classname = classes{i};
                nrobj = length(this_app.objects.(classname).data);
                for k = 1 : nrobj
                    h = plotobject(this_app,this_app.objects.(classname).data{k},this_app.objects.(classname).color);
                    this_app.objects.(classname).handles(k) = h;
                    if ~this_app.objects.(classname).visible(k)
                        set(h,'Visible','off');
                    end
                end
            end
            hold off
        else
            
            % create new instance of topoapp
            validateattributes(this_app.DEM,{'GRIDobj'},{},'topoapp','DEM',1)

            tic
            fprintf(1,'Initializing topoapp...');
            % create coordinate vectors
            [this_app.X this_app.Y] =  getcoordinates(this_app.DEM);
            this_app.X = this_app.X(:)';
            this_app.Y = this_app.Y(:);
            % create slope map
            this_app.G = gradient8(this_app.DEM,'degree');
            % create elevation-colored hillshade image
            this_app.RGB = imageschs(this_app.DEM,this_app.DEM);
            % create grayscale hillshade image
            this_app.HS = imageschs(this_app.DEM,nan(this_app.DEM.size));
            % create slope-colored hillshade image
            this_app.HSG = imageschs(this_app.DEM,this_app.G,'colormap','jet');
            % create hillshade with sink depth image
            SINKDEPTH = fillsinks(this_app.DEM)-this_app.DEM;
            SINKDEPTH.Z(SINKDEPTH.Z==0) = nan;
            this_app.SINKS   = imageschs(this_app.DEM,SINKDEPTH); 
            this_app.objects.classes = '';
            fprintf(1,'done. ');
            toc
            
            guilayout(this_app) % create graphical user interface
            this_app = defaultparas(this_app);
            
            
        end
        
    end % 
    %----------------------------------------------------------------
    function guilayout(app)
        
        tic
        fprintf(1,'Initializing gui...');
        % Create DEM window
        scrsz = get(0,'ScreenSize'); % = [left, bottom, width, height] in pixels
        hgt = scrsz(4)*2/3; ratio = app.DEM.size(2)/app.DEM.size(1)-0.1;
        app.gui.hfig = figure('Name','DEMfigure','OuterPosition',[1/4*scrsz(3) 1/3*scrsz(4) hgt*ratio  hgt]);
        app.gui.hax = imgca(app.gui.hfig);
        app.gui.himg = image(app.RGB,'parent',app.gui.hax,'Tag','DEM');
        
        % Setup interactive zooming
        set(app.gui.hax,'ylimmode','auto','xlimmode','auto'); % reset after image detroyed
        app.gui.hzoom = zoom(app.gui.hax); % build zoom object
        zoom reset; % store current setting
        
        % Create an instance of mscrollpanel and use API
        app.gui.hscroll = imscrollpanel(app.gui.hfig,app.gui.himg);
        api = iptgetapi(app.gui.hscroll);
        mag = api.findFitMag();
        api.setMagnification(mag); % set to initial magnification
        
        % Magnification window
        app.gui.hmagbox = immagbox(app.gui.hfig,app.gui.himg);
        pos = get(app.gui.hmagbox,'Position');
        set(app.gui.hmagbox,'Position',[0 0 pos(3) pos(4)])
        imoverview(app.gui.himg)
        fprintf(1,'done. ');
        toc
            
        % Setup toolbar
        app.gui.hTB = uitoolbar(app.gui.hfig);
        
        % Initilize tools
        tic
        fprintf(1,'Initializing tools:\n');
        app.gui.TB = [];
        sep = '\'; if ismac; sep = '/'; end
        A = what('@topoapp');
        fid = fopen([A.path,sep,'topoapptools.txt']);
        C = textscan(fid,'%s'); fclose(fid);
        C = deblank(C{1});
        ct=1; addsep=0;
        while ct<=length(C)
            if ~isempty(C{ct})
                if strcmp(C{ct},'---')
                    addsep=1;
                else
                    disp(C{ct})
                    [app] = feval(C{ct},0,'init',app);
                    if addsep
                        set(app.gui.TB(end),'Separator','on')
                        addsep=0;
                    end
                end
            end
            ct=ct+1;
        end
        toc
        
    end % 
    %----------------------------------------------------------------
    function app = defaultparas(app)
        % flowrouting parameters
        app.paras.flowset.minarea = 2e6;  % minimum upstream area for drainage network [km^2]
        app.paras.flowset.shorties = 1e3; % maximum length of first order streams to remove [m]
    end
    %----------------------------------------------------------------
    
    
end % methods

end % classdef

