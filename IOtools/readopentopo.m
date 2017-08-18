function DEM = readopentopo(varargin)

%READOPENTOPO read DEM using the opentopography.org API
%
% Syntax
%
%     DEM = readopentopo(pn,pv,...)
%
% Description
%
%     readopentopo reads DEMs from opentopography.org using the API
%     described on:
%     http://www.opentopography.org/developers
%     The DEM comes in geographic coordinates (WGS84) and should be
%     projected to a projected coordinate system (use reproject2utm) before
%     analysis in TopoToolbox.
%
% Input arguments
%
%     Parameter name values
%     'filename'       provide filename. By default, the function will save
%                      the DEM to a temporary file in the system's temporary 
%                      folder.
%     'north'          northern boundary in geographic coordinates (WGS84)
%     'south'          southern boundary
%     'west'           western boundary
%     'east'           eastern boundary
%     'demtype'        The global raster dataset - SRTM GL3 (90m) is 
%                      'SRTMGL3', SRTM GL1 (30m) is 'SRTMGL1', SRTM GL1 
%                      (Ellipsoidal) is 'SRTMGL1_E', and ALOS World 3D 30m  
%                      is 'AW3D30'
%     'deletefile'     'true' or false. True, if file should be deleted
%                      after it was downloaded and added to the workspace.
% 
% Output arguments
%
%     DEM            Digital elevation model in geographic coordinates
%                    (GRIDobj)
%
% See also: GRIDobj, websave
%
% Reference: http://www.opentopography.org/developers
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. June, 2017


p = inputParser;
addParameter(p,'filename',[tempname '.tif']);
% addParameter(p,'interactive',false);
addParameter(p,'north',37.091337);
addParameter(p,'south',36.738884);
addParameter(p,'west',-120.168457);
addParameter(p,'east',-119.465576);
addParameter(p,'demtype','SRTMGL3');
addParameter(p,'deletefile',true);
parse(p,varargin{:});

demtype = validatestring(p.Results.demtype,{'SRTMGL3','SRTMGL1','SRTMGL1_E','AW3D30'},'readopentopo');

url = 'http://opentopo.sdsc.edu/otr/getdem';

% create output file
f = fullfile(p.Results.filename);

% save to drive
options = weboptions('Timeout',inf);

west = p.Results.west;
east = p.Results.east;
south = p.Results.south;
north = p.Results.north;

% if any([isempty(west) isempty(east) isempty(south) isempty(north)]) || p.Results.interactive;
%     wm = webmap;
%     % get dialog box
%     messagetext = ['Zoom and resize the webmap window to choose DEM extent. ' ...
%                          'Click the close button when you''re done.'];
%     d = waitdialog(messagetext);
%     uiwait(d);    
%     [latlim,lonlim] = wmlimits(wm);
%     west = lonlim(1);
%     east = lonlim(2);
%     south = latlim(1);
%     north = latlim(2);
% end
    
    

websave(f,url,'west',west,...
              'east',east,...
              'north',north,...
              'south',south,...
              'outputFormat', 'GTiff', ...
              'demtype', demtype, ...
              options);

DEM      = GRIDobj(f);
DEM.name = demtype;

if p.Results.deletefile
    delete(f);
end
end

% function d = waitdialog(messagetext)
%     d = dialog('Position',[300 300 250 150],'Name','Choose rectangle region',...
%         'WindowStyle','normal');
% 
%     txt = uicontrol('Parent',d,...
%                'Style','text',...
%                'Position',[20 80 210 40],...
%                'String',messagetext);
% 
%     btn = uicontrol('Parent',d,...
%                'Position',[85 20 70 25],...
%                'String','Close',...
%                'Callback','delete(gcf)');
% end


