function [GS,val] = STREAMobj2shape(S,varargin)

%STREAMobj2shape Convert STREAMobj to geoshape or mapshape
%
% Syntax
%
%     GS = STREAMobj2shape(S)
%     GS = STREAMobj2shape(S,'seglength',seglength,'type',type)
%     GS = STREAMobj2shape(...,'attributes',{attname data aggfun ...})
%     GS = STREAMobj2shape(...,'summarizeby',attname)
%
% Description
%
%     STREAMobj2shape takes a STREAMobj and converts it to a geoshape or
%     mapshape object. geoshape or mapshape objects are ways to store
%     vector features with either point, line, or polygon topology. They
%     can be exported as shapefile, and hopefully soon support projections.
%
% Input parameters
%
%     S       STREAMobj
%     
%     Parameter name/value pairs
%
%     'seglength'   length of line features
%     'type'        'geo' or 'map'. Determines whether the function returns
%                   a geoshape or mapshape object.
%     'attributes'  cell array with attribute data (see STREAMobj2mapstruct
%                   for details). These data will be aggregated to feature
%                   properties. 
%
% Output arguments
%
%     GS      geoshape or mapshape object
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     c = chitransform(S,A);
% 
%     GS = STREAMobj2shape(S,'type','geo',...
%             'attributes',{...
%             'z' DEM @mean ...
%             'chi' c @mean});
%     colorRange = makesymbolspec('Line',...
%             {'chi',[min(c) max(c)],'Color',jet(100)});
%     geoshow(GS,'symbolspec',colorRange)
%
% Example 2 Summarize by attribute
%
%     [GS,val] = STREAMobj2shape(S,'type','geo',...
%             'attributes',{...
%               'z' DEM @mean ...
%               'chi' c @mean},...
%             'summarizeby','chi');
%     clr = jet(numel(val));
%     for r = 1:numel(val); 
%          geoshow(GS{r},'color',clr(r,:));
%     end
%
% Example 3 Plot on webmap
%
%     for r = 1:numel(val); 
%          wmline(GS{r},'color',clr(r,:));
%     end
% 
% See also: STREAMobj2mapstruct
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 16. October, 2019

d = S.distance;
defaultseglength = max(max(d)/20,5*S.cellsize);

p = inputParser;
addParameter(p,'type','mapshape',...
    @(x) ischar(validatestring(x,{'geoshape' 'mapshape'},'STREAMobj2shape')));
addParameter(p,'seglength',defaultseglength,@(x) isscalar(x) && x>S.cellsize);
addParameter(p,'attributes',{})
addParameter(p,'summarizeby','')
addParameter(p,'by',20)
parse(p,varargin{:});

MS = STREAMobj2mapstruct(S,'seglength',p.Results.seglength,...
    'attributes',p.Results.attributes);

% remove nans from coordinates if any


if ~isempty(p.Results.summarizeby)
    if ~isfield(MS,p.Results.summarizeby)
        error('Invalid summary field')
    end
end

type = validatestring(p.Results.type,{'geoshape' 'mapshape'});

switch type
    case 'geoshape'
        [lat,lon] = cellfun(@(x,y) minvtran(S.georef.mstruct,x,y),{MS.X},{MS.Y},'UniformOutput',false);
        [MS.Lat] = deal(lat{:});
        [MS.Lon] = deal(lon{:});
        
        MS = rmfield(MS,{'X' 'Y'});
end

if ~isempty(p.Results.summarizeby)
    [~,edges,bins] = histcounts([MS.(p.Results.summarizeby)],p.Results.by);
    split = true;
else
    split = false;
end


switch type
    case 'geoshape'
        if split
            GS = accumarray(bins(:),(1:numel(MS))',[],@(ix) {geoshape(MS(ix))});
            val = edges(1:end-1) + diff(edges)/2;
        else
            GS = geoshape(MS);
            val = [];
        end
    case 'mapshape'
        if split
            GS = accumarray(bins(:),(1:numel(MS))',[],@(ix) geoshape(MS(ix)));
            val = edges(1:end-1) + diff(edges)/2;
        else
            GS = mapshape(MS);
            val = [];
            try 
                GS.GTCitationGeoKey = DEM.georef.GeoKeyDirectoryTag.GTCitationGeoKey;
            end
        end  
end
    
    
