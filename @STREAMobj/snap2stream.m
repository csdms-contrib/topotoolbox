function [xn,yn,IX,D,MP] = snap2stream(S,x,y,varargin)

%SNAP2STREAM snap locations to nearest stream location
%
% Syntax
%
%     [xn,yn] = snap2stream(S,x,y)
%     [xn,yn] = snap2stream(S,x,y,pn,pv,...)
%     [xn,yn] = snap2stream(S,lat,lon,'inputislatlon',true,pn,pv,...)
%     [xn,yn,IX,res,MP] = snap2stream(...)
%
% Description
%
%     Often gauge locations are not precisely associated with 
%     streams derived from a digital elevation model. Calculation of 
%     drainage basins may thus return wrong basin delineations. One method 
%     is to manually adjust their location to a derived stream raster.
%     snap2stream does this automatically by finding identifying the
%     shortest distance to a stream. 
%
% Input
% 
%     S      instance of STREAMobj
%     x,y    coordinate vectors
%     
%     Parameter name/value pairs {default}
%
%     'snapto'   string 
%     constrain snapping to network features. Features can be {'all'},
%     'confluence', 'outlet' or 'channelhead'
%
%     'streamorder'   scalar or string
%     scalar: if you supply a scalar integer s as parameter value,
%     locations will be snapped only to streams with order s. string: a
%     string allows you to define a simple relational statement, e.g. '>=3'
%     will snap only to streams of streamorder greater or equal than 3.
%     Other commands may be '==5', '<5', etc.
%
%     'maxdist'  scalar
%     maximum snapping distance {inf}. 
%
%     'inputislatlon'   true or {false}
%     determines whether input coordinates are supplied in a geographic
%     coordinate system (wgs84). If set to true, snap2stream projects the
%     geographic coordinates to the same coordinate system as S. This
%     requires that S has a projected coordinate system and the mapping
%     toolbox. 
%
%     'plot'  true or {false}
%     plot results
%
%
% Output arguments
%
%     xn,yn    coordinate vectors of snapped locations
%     IX       linear index into grid from which S was derived
%     res      residual (eucl. distance between locations)
%     MP       map struct that can be exported to a shapefile
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. January, 2018


% Parse Inputs
p = inputParser;         
p.FunctionName = 'snap2stream';

snaptomethods = {'all','confluence','outlet','channelhead'};

addParamValue(p,'snapto','all',@(x) ischar(validatestring(x,snaptomethods)));
addParamValue(p,'maxdist',inf,@(x) isscalar(x) && x>0);
addParamValue(p,'streamorder',[]);
addParamValue(p,'inputislatlon',false,@(x) isscalar(x));
addParamValue(p,'plot',false,@(x) isscalar(x));

parse(p,varargin{:});

snapto = validatestring(p.Results.snapto,{'all','confluence','outlet','channelhead'});
I = true(size(S.x));

if p.Results.inputislatlon
    if isempty(S.georef)
        error('TopoToolbox:georeferencing',...
            ['Projection of stream network unknown. Cannot convert ' ...
             'geographic coordinates.']);
    end
    [x,y] = mfwdtran(S.georef.mstruct,x,y);
end

% evaluate optional arguments to restrict search to specific stream
% locations only

switch snapto;
    case 'all'
        % do nothing
    otherwise
        nrc = numel(S.x);
        M = sparse(S.ix,S.ixc,true,nrc,nrc);
        
        switch snapto;
            case 'confluence'
                I = I & (sum(M,1)>1)';
            case 'outlet'
                I = I & (sum(M,2) == 0);
            case 'channelhead'
                I = I & (sum(M,1) == 0)';
        end
end
        
if ~isempty(p.Results.streamorder);
    
    so = streamorder(S);
    if isnumeric(p.Results.streamorder)
        validateattributes(p.Results.streamorder,{'numeric'},...
            {'scalar','integer','>',0,'<=',max(so)});
        I = (so == p.Results.streamorder) & I;
    else
        try
            % expecting relational operator
            sothres = nan;
            counter = 1;
            while isnan(sothres);
                sothres = str2double(p.Results.streamorder(counter+1:end));
                relop   = p.Results.streamorder(1:counter);
                counter = counter+1;
            end
        catch ME
            rethrow(ME)
        end

        relop = str2func(relop);
        I = relop(so,sothres) & I;
        
    end
end

% find nearest stream location

xn = S.x(I);
yn = S.y(I);
IX = S.IXgrid(I);

xn = xn(:);
yn = yn(:);
IX = IX(:);

x  = x(:);
y  = y(:);

% fast and easy using Statistics toolbox
try
    [IDX,D] = knnsearch([xn yn],[x y]);
catch %#ok<CTCH>
    % aargh, not available... dirty implementation
    D = nan(numel(x),1);
    IDX = nan(numel(x),1);
    
    for r = 1:numel(x);
        [D(r), IDX(r)] = min(sqrt((xn-x(r)).^2 + (yn-y(r)).^2));
    end           
end

xn = xn(IDX);
yn = yn(IDX);
IX = IX(IDX);



% apply maximum distance
if ~isinf(p.Results.maxdist);
    ID = D <= p.Results.maxdist;
else
    ID = true(numel(x),1);
end

if nargout == 5;
    % create output mapstruct
    MP = struct('Geometry',{'Point'},...
        'X',num2cell(xn),...
        'Y',num2cell(yn),...
        'IX',num2cell(IX),...
        'residual',num2cell(D),...
        'inrange',num2cell(double(ID)),...
        'Xold',num2cell(x),...
        'Yold',num2cell(y));
end



if ~isinf(p.Results.maxdist);
    ID = ~ID;
    xn(ID) = nan;
    yn(ID) = nan;
    IX(ID) = nan;
    ID = ~ID;
end

% plot results
if p.Results.plot
    plot(S)
    hold on
    plot(x,y,'sk')
    plot(xn(ID),yn(ID),'xr')
    plot([xn(ID) x(ID)]',[yn(ID) y(ID)]','k');
end
    




