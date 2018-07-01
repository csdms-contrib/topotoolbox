function [cmap,zlimits] = ttcmap(zlimits,varargin)

%TTCMAP create elevation color map optimized for elevation range
%
% Syntax
%     
%     [cmap,zlimits] = ttcmap(zlimits)
%     [cmap,zlimits] = ttcmap(DEM)
%     [cmap,zlimits] = ttcmap(...,pn,pv,...)
%     ttcmap
%     t = ttcmap
%
% Description
%
%     TTCMAP has a functionality similar to the Mapping Toolbox function
%     demcmap. The function creates a colormap for displaying digital
%     elevation models, in particularly if topography and bathymetry are to
%     be displayed in the same map. ttcmap returns a colormap so that
%     negative values in the DEM are distinct from positive values. To
%     ensure that the color transition is exactly at the value 0, the color
%     limits (e.g. via the command clim) in the displayed map must be set to 
%     zlimits. 
%
%     ttcmap without in- and output arguments shows available colormaps and
%     their elevation range.
%
%     t = ttcmap with one output argument but no input arguments returns
%     available colormaps and their elevation range as table.
%
% Input arguments
%
%     zlimits        two element vector with maximum and minimum elevation
%     DEM            GRIDobj from which zlimits will be calculated
%    
%     Parameter name/value pairs
%
%     'n'       number of colors (default = 255).
%     'cmap'    colormap to be used. Default is 'gmtrelief'. For an overview
%               of available colormaps call ttcmap with no input arguments
%     'zero'    'land' or 'sea'. Determines whether the value zero will be
%               displayed as 'land' (default) or as 'sea'.
%     
%     
% Output arguments
%
%     cmap      n*3 colormap
%     zlimits   adjusted elevation range 
%     t         table with available colormaps and their elevation ranges
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     % Let's pretend that the lower 300 m of the DEM
%     % are submerged
%     DEM = DEM-min(DEM)-300;
% 
%     [clr,zlimits] = ttcmap(DEM,'cmap','gmtglobe');
%     imageschs(DEM,DEM,'colormap',clr,'caxis',zlimits);
%
%     % If you want to plot with imagesc
%     figure
%     imagesc(DEM);
%     colormap(clr)
%     caxis(zlimits)
%
%
% References: Colormaps available through ttcmap are obtained from 
%     following website:
%     http://soliton.vm.bytemark.co.uk/pub/cpt-city/views/totp-svg.html
%
% See also: IMAGESCHS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. February, 2018

allowedcmaps = {'gmtrelief','france','mby','gmtglobe','etopo1'};
if nargin == 0
    ncmaps = numel(allowedcmaps);
    Min_Elevation = zeros(ncmaps,1);
    Max_Elevation = zeros(ncmaps,1);
    for r = 1:numel(allowedcmaps)
        elev = getcolormap(allowedcmaps{r}); 
        Min_Elevation(r) = min(elev);
        Max_Elevation(r) = max(elev);
    end

    Colormap_Name = allowedcmaps(:);
    cmap = table(Colormap_Name,Min_Elevation,Max_Elevation);
    
    if nargout == 2
        zlimits = [];
    elseif nargout == 0
        disp(cmap);
        clear cmap
        return
    end
           
    return
end

p = inputParser;
p.FunctionName = 'ttcmap';
% optional
addParamValue(p,'n',255,@(x)isnumeric(x) && x>=1);
addParamValue(p,'cmap','gmtrelief');
addParamValue(p,'zero','land');
parse(p,varargin{:});

cmaptype = validatestring(p.Results.cmap,allowedcmaps);

if isa(zlimits,'GRIDobj')
    zlimits = [min(zlimits) max(zlimits)];
end

nr = 255;
[elev,clr] = getcolormap(cmaptype);


% Problem: the color map should 
elevmap = linspace(min(zlimits(:)),max(zlimits(:)),nr);
% which of the entries is closest to zero
if zlimits(1)<0 && zlimits(2)>0
    ix = find(elev==0);
    
    if ~isempty(ix)
        
        switch lower(p.Results.zero)
            case 'land'
                elev = [elev(1:ix-1);elev(ix-1)./1000; elev(ix:end)];
                clr  = [clr(1:ix-1,:);clr(ix-1,:); clr(ix:end,:)];
            case 'water'
                elev = [elev(1:ix);elev(ix)./1000; elev(ix+1:end)];
                clr  = [clr(1:ix,:);clr(ix,:); clr(ix+1:end,:)];
        end
    end
            
    [~,ix] = min(abs(elevmap));
    elevmap = elevmap-elevmap(ix);
    zlimits = [elevmap(1) elevmap(end)];
else
    zlimits = [min(zlimits) max(zlimits)];
end

cmap = interp1(elev,clr,elevmap);

% if the elevation range in zlimits extends beyond the limits of the
% colormap, interp1 will return nans in the colormap
I = isnan(cmap);
if all(I(:,1))
    error('The zlimits are outside the range of the colormap')
end
if any(I(:,1))
    I = ~I(:,1);
    ix = find(I,1,'first');
    for r = 1:3; cmap(1:ix,r) = cmap(ix,r);end
    ix = find(I,1,'last');
    for r = 1:3; cmap(ix:end,r) = cmap(ix,r);end
end
        


end

function [elev,clr] = getcolormap(type)

switch lower(type)
    
    %% gmtglobe
    case 'gmtglobe'
    
data = ...
[-10000	153	0	255
-9500	153	0	255
-9000	153	0	255
-8500	136	17	255  
-8000	119	34	255
-7500	102	51	255
-7000	85	68	255
-6500	68	85	255
-6000	51	102	255
-5500	34	119	255  
-5000	17	136	255
-4500	0	153	255
-4000	27	164	255
-3500	54	175	255
-3000	81	186	255
-2500	108	197	255
-2000	134	208	255
-1500	161	219	255
-1000	188	230	255
-500	215	241	255
-200	241	252	255
0	51	102	0
100	51	204	102
200	187	228	146
500	255	220	185
1000	243	202	137
1500	230	184	88
2000	217	166	39
2500	168	154	31
3000	164	144	25
3500	162	134	19
4000	159	123	13
4500	156	113	7
5000	153	102	0
5500	162	89	89
6000	178	118	118
6500	183	147	147
7000	194	176	176
7500	204	204	204
8000	229	229	229
8500	242	242	242
9000	255	255	255
9500	255	255	255];
    case 'france'
        %% france
data = [...
-5000 113 171 216
-4500 121 178 222
-3500 132 185 227
-2500 141 193 234
-2000 150 201 240
-1500 161 210 247
-1000 172 219 251
-500  185 227 255
-200  198 236 255
-100  216 242 254
0     172 208 165
100   148 191 139
250   189 204 150
500   225 228 181
750   239 235 192
1000  222 214 163
1500  202 185 130
2000  172 154 124
2500  186 174 154
3000  202 195 184
3500  224 222 216
4000  245 244 242];
    case 'gmtrelief'
        data = [...
-8000	0	0	0	
-7000	0	5	25	
-6000	0	10	50	
-5000	0	80	125	
-4000	0	150	200	
-3000	86	197	184	
-2000	172	245	168	
-1000	211	250	211	
0	70	120	50
500	120	100	50
1000	146	126	60	
2000	198	178	80	
3000	250	230	100	
4000	250	234	126	
5000	252	238	152	
6000	252	243	177	
7000	253	249	216	];
    case 'mby'
data = [...
-8000   0   0	80  
-6000   0   30  100 
-4000   0   50  102 
-2500	19  108 160	
-150	24  140 205	
-50     135 206 250 
-10	    176 226	255	
0       0   97	71  
50	16  123	48
500	232 214	125	
1200	163 68  0	
1800	130 30  30	
2800	161 161 161	
4000	206 206 206	];

    case 'etopo1'
data = ...
[-11000	10	0	121
-10500	26	0	137
-10000	38	0	152
-9500	27	3	166
-9000	16	6	180
-8500	5	9	193
-8000	0	14	203
-7500	0	22	210
-7000	0	30	216
-6500	0	39	223
-6000	12	68	231
-5500	26	102	240
-5000	19	117	244
-4500	14	133	249
-4000	21	158	252
-3500	30	178	255
-3000	43	186	255
-2500	55	193	255
-2000	65	200	255
-1500	79	210	255
-1000	94	223	255
-500	138	227	255
0.00	51	102	0
100	    51	204	102
200	    187	228	146
500	    255	220	185
1000	243	202	137
1500	230	184	88
2000	217	166	39
2500	168	154	31
3000	164	144	25
3500	162	134	19
4000	159	123	13
4500	156	113	7
5000	153	102	0
5500	162	89	89
6000	178	118	118
6500	183	147	147
7000	194	176	176
7500	204	204	204
8000	229	229	229];

    otherwise
        error('unknown colormap');
end

elev = data(:,1);
clr  = data(:,2:end)/255;
end
