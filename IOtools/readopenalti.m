function data = readopenalti(varargin)

%READOPENALTI read altimetry data using the openaltimetry.org API
%
% Syntax
%
%     data = readopenalti(pn,pv,...)
%
% Description
%
%     readopenalti reads altimetry data from openaltimetry.org using the 
%     API described on: https://openaltimetry.org/data/swagger-ui/
%     The data comes in geographic coordinates (WGS84). EGM96 geoid heights
%     are added using the Mapping Toolbox function egm96geoid 
%
% Input arguments
%
%     Parameter name values
%     'filename'       provide filename. By default, the function will save
%                      the data to a temporary file in the system's temporary 
%                      folder. The option 'deletefile' controls whether the
%                      file is kept on the hard drive.
%     'date'           datetime scalar or vector indicating the days
%                      required. Default is the last 600 days.
%     'extent'         GRIDobj or four element vector with geographical 
%                      coordinates in the order [west east south north].
%                      If a GRIDobj is supplied, readopentopo uses the
%                      function GRIDobj/getextent to obtain the bounding
%                      box in geographical coordinates. If extent is set,
%                      then the following parameter names 'north',
%                      'south', ... are ignored.
%     'addmargin'      Expand the extent derived from 'extent',GRIDobj by a
%                      scalar value in °. Default is 0.01. The option is
%                      only applicable if extent is provided by a GRIDobj.
%     'north'          northern boundary in geographic coordinates (WGS84)
%     'south'          southern boundary
%     'west'           western boundary
%     'east'           eastern boundary
%     'product'        'atl03' L2A Global Geolocated Photon Data
%                      'atl06' L3A Land Ice Height
%                      'atl07' L3A Sea Ice Height
%                      'atl08' L3A Land and Vegetation Height (default)
%                      'atl10' L3A Sea Ice Freeboard
%                      'atl12' L3A Ocean Surface Height
%                      'atl13' L3A Inland Water Surface Height, Version 1                     
%     'level3a'        {true} or false
%     'verbose'        {true} or false. If true, then some information on
%                      the process is shown in the command window
%     'deletefile'     {true} or false. True, if file should be deleted
%                      after it was downloaded and added to the workspace.
% 
% Output arguments
%
%     data     table with altimetry data
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     d = readopenalti('extent',DEM,'addmargin',0,...
%           'date',datetime(2020,05,1):datetime('today'));
%     % coordinates are returned as longitude and latitude
%     [d.x,d.y] = mfwdtran(DEM.georef.mstruct,d.latitude,d.longitude);
%     surf(DEM,'block',true);
%     colormap(landcolor)
%     camlight
%     axis off
%     hold on
%     % plot points with some offset to elevation so that points are
%     % plotted above the surface
%     plot3(d.x,d.y,d.h_te_best_fit+20,'.r')
%
% See also: GRIDobj, websave, readopentopo, egm96geoid
%
% Reference: https://openaltimetry.org/data/swagger-ui/
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. November, 2021


defaultdate = datetime('today')-600 : datetime('today');

p = inputParser;
addParameter(p,'filename',[tempname '.csv']);
addParameter(p,'interactive',false);
addParameter(p,'extent',[]);
addParameter(p,'addmargin',0.01);
addParameter(p,'north',37.091337);
addParameter(p,'south',36.738884);
addParameter(p,'west',-120.168457);
addParameter(p,'east',-119.465576);
addParameter(p,'date',defaultdate)
addParameter(p,'product','atl08');
addParameter(p,'level3a',false);
addParameter(p,'deletefile',true);
addParameter(p,'verbose',true);
parse(p,varargin{:});

validproducts = {'atl03','atl10','atl12',...
                 'atl13','atl06','atl07',...
                 'atl08'};

product = validatestring(p.Results.product,validproducts,'readopenalti');
urltracks = 'https://openaltimetry.org/data/api/icesat2/getTracks?';
if ~p.Results.level3a
    url = ['https://openaltimetry.org/data/api/icesat2/' product '?'];
else
    url = 'https://openaltimetry.org/data/api/icesat2/level3a?';
end
% Example call
% https://openaltimetry.org/data/api/icesat2/atl08?date=2020-05-13&minx=55&miny=46&maxx=57&maxy=48&trackId=730&client=portal&outputFormat=csv


% create output file
f = fullfile(p.Results.filename);
    

% save to drive
options = weboptions('Timeout',100000);

% get extent
if ~isempty(p.Results.extent)
    if isa(p.Results.extent,'GRIDobj')
        ext = getextent(p.Results.extent,true);
        west  = ext(1) - p.Results.addmargin;
        east  = ext(2) + p.Results.addmargin;
        south = ext(3) - p.Results.addmargin;
        north = ext(4) + p.Results.addmargin;
        
    elseif numel(p.Results.extent) == 4
        west = p.Results.extent(1);
        east = p.Results.extent(2);
        south = p.Results.extent(3);
        north = p.Results.extent(4);
    else
        error('Unknown format of extent')
    end
else
    west = p.Results.west;
    east = p.Results.east;
    south = p.Results.south;
    north = p.Results.north;
end

% now we have an extent. Or did the user request interactively choosing
% the extent.
if any([isempty(west) isempty(east) isempty(south) isempty(north)]) || p.Results.interactive
    wm = webmap;
    % get dialog box
    messagetext = ['Zoom and resize the webmap window to choose DEM extent. ' ...
                         'Click the close button when you''re done.'];
    d = waitdialog(messagetext);
    uiwait(d);    
    [latlim,lonlim] = wmlimits(wm);
    west = lonlim(1);
    east = lonlim(2);
    south = latlim(1);
    north = latlim(2);
end
    
if p.Results.verbose
    a = areaint([south south north north],...
                [west east east west],almanac('earth','radius','kilometers'));
    disp('-------------------------------------')
    disp('readopenalti process:')
    disp(['Product product: ' product])
    disp(['API url tracks: ' urltracks])
    disp(['API url product: ' url])
    disp(['Local file name: ' f])
    disp(['Area: ' num2str(a,2) ' sqkm'])
    disp('-------------------------------------')
    disp(['Starting download: ' datestr(now)])
end

% get date
if isdatetime(p.Results.date)
    dt = p.Results.date;
    D  = string(dt,'yyyy-MM-dd');
else
    D = p.Results.date;
    dt = datetime(D,'InputFormat','yyyy-MM-dd');
end

if p.Results.verbose
    disp(['Get tracks: ' datestr(now)])
end

% Download tracks with websave
tracktable = table([],[],'VariableNames',{'track','date'});
if p.Results.verbose
    disp(['Identifying tracks: ' datestr(now)])
end

for r=1:numel(dt)
	   textwaitbar(r,numel(dt),'Wait');
	   
	   outfile = websave(f,urltracks,...
              'date',D(r),...
              'minx',west,...
              'maxx',east,...
              'maxy',north,...
              'miny',south,...
              'client','TopoToolbox',...
              'outputFormat', 'csv', ...
              options);
    temptracks = readtable(f);
    temptracks.date = repmat(dt(r),size(temptracks,1),1);
    tracktable = [tracktable; temptracks];
    delete(f);
	   
end


for r = 1:numel(dt)
    
end

totaltracksfound = size(tracktable,1);

if p.Results.verbose
    disp([num2str(totaltracksfound) ' tracks found.'])
    disp(['Data download starts: ' datestr(now)])
end

% Download data
if p.Results.level3a
    
    counter = 1;
    trackIds = unique(tracktable.track);
    startDate = D(1);
    endDate = D(end);
    for r = 1:numel(trackIds)
        trackId = trackIds(r);
        if p.Results.verbose
            disp(['Download trackId ' num2str(trackId) ': ' datestr(now)])
        end
        try
            outfile = websave(f,url,...
                'product',product,...
                'startDate',startDate,...
                'endDate',endDate,...
                'minx',west,...
                'maxx',east,...
                'maxy',north,...
                'miny',south,...
                'trackId',trackId,...
                'client','TopoToolbox',...
                'outputFormat', 'csv', ...
                options);
            if counter == 1 || ~exist('data','var')
                data = readtable(f);
            else
                data = [data; readtable(f)];
            end
            counter = counter + 1;
        catch
            disp(['Download failed: ' datestr(now)])
        end
        delete(f);
    end
    
else
    
    counter = 1;
    for r = 1:totaltracksfound
        dd = string(tracktable.date(r),'yyyy-MM-dd');
        trackId = num2str(tracktable.track(r));
        if p.Results.verbose
            disp(['Download ' char(dd) ',' num2str(trackId) ': ' datestr(now)])
        end
        try
            outfile = websave(f,url,...
                'date',dd,...
                'minx',west,...
                'maxx',east,...
                'maxy',north,...
                'miny',south,...
                'trackId',trackId,...
                'client','TopoToolbox',...
                'outputFormat', 'csv', ...
                options);
            if counter == 1  || ~exist('data','var')
                data = readtable(f);
            else
                data = [data; readtable(f)];
            end
            counter = counter + 1;
        catch
            disp(['Download failed: ' datestr(now)])
        end
        delete(f);
    end
end

if p.Results.verbose
    disp(['Total of ' num2str(size(data,1)) ' points downloaded: ' datestr(now)])

end

if p.Results.verbose
    disp(['Adding EGM96 heights: ' datestr(now)]);
end

data.egm96geoid = egm96geoid(data.latitude,data.longitude);

if p.Results.verbose
    disp(['Done: ' datestr(now)])
    disp('-------------------------------------')
end
end

function d = waitdialog(messagetext)
    d = dialog('Position',[300 300 250 150],'Name','Choose rectangle region',...
        'WindowStyle','normal');

    txt = uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 210 40],...
               'String',messagetext);

    btn = uicontrol('Parent',d,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
end

function textwaitbar(i, n, msg)
% A command line version of waitbar.
% Usage:
%   textwaitbar(i, n, msg)
% Input:
%   i   :   i-th iteration.
%   n   :   total iterations.
%   msg :   text message to print.
%
% Date      : 05/23/2019
% Author    : Xiaoxuan He   <hexxx937@umn.edu>
% Institute : University of Minnesota
%
% Previous percentage number.
persistent i_prev_prct;
% Current percentage number.
i_prct = floor(i ./ n * 100);
% Print message when counting starts.
if isempty(i_prev_prct) || i_prct < i_prev_prct
    i_prev_prct = 0;
    S_prev = getPrctStr(i_prev_prct);
    
    fprintf('%s: %s',msg, S_prev);
end
% Print updated percentage.
if i_prct ~= i_prev_prct
    S_prev = getPrctStr(i_prev_prct);
    fprintf(getBackspaceStr(numel(S_prev)));
    
    S = getPrctStr(i_prct);
    fprintf('%s', S);
    
    i_prev_prct = i_prct;
end
% Clear percentage variable.
if i_prct == 100
    fprintf(' Done.\n');
    clear i_prev_prct;
end
end
function S = getPrctStr(prct)
S = sprintf('%d%%  %s',prct,getDotStr(prct));
if prct < 10
    S = ['  ',S];
elseif prct < 100
    S = [' ',S];
end
end
function S = getDotStr(prct)
S = repmat(' ',1,10);
S(1:floor(prct/10)) = '.';
S = ['[',S,']'];
end
function S = getBackspaceStr(N)
S = repmat('\b',1,N);
end

