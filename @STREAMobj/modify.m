function S = modify(S,varargin)

% modify instance of STREAMobj to meet user-defined criteria
%
% Syntax
%
%     S2 = modify(S,pn,pv)
%
% Description
%
%     The function modify changes the geometry of a stream network to meet
%     different user-defined criteria. See
%
%     demo_modifystreamnet
%
%     for an overview of the function's scope.
%
% Input arguments
%
%     S     instance of stream object
%
%     Parameter name/value pairs (one pn/pv pair is required to run the
%     function)
%
%     'streamorder'   scalar or string
%     scalar: if you supply a scalar integer x as parameter value, S2
%     contains only streams with order x.
%     string: a string allows you to define a simple relational statement,
%     e.g. '>=3' returns the stream network with streams of streamorder
%     greater or equal than 3. Other commands may be '==5', '<5', etc.
%
%     'distance' scalar or [2x1] vector
%     lets you define a minimum and maximum distance from which the stream 
%     network starts and ends based on the distances given in S. A scalar is
%     interpreted as minimum upstream distance and a two element vector is
%     interpreted as minimum and maximum upstream distance to which the
%     stream network is cut.
%
%     'maxdsdistance' scalar
%     modifies the stream network that only the portion of the network is
%     retained that is within downstream distance of the the channelheads.
%
%     'upstreamto' logical GRIDobj or linear index in GRIDobj
%     returns the stream network upstream to true pixels in the logical
%     raster of GRIDobj. Note that, if the grid contains linear features
%     (e.g. a fault), the line should be 4 connected. Use
%     bwmorph(I.Z,'diag') to establish 4-connectivity in a logical raster
%     I.
%
%     'downstreamto' logical GRIDobj or linear index in GRIDobj
%     returns the stream network downstream to true pixels in the logical
%     raster of GRIDobj. Note that, if the grid contains linear features
%     (e.g. a fault), the foreground (true pixels) should be 4 connected.
%     Use bwmorph(I.Z,'diag') to establish 4-connectivity in a logical
%     raster I.
%
%     'rmconncomps' scalar
%     removes connected components (individual stream 'trees') of the
%     entire network with a maximum distance in map units less than the 
%     specified value.
%
%     'rmconncomps_ch' scalar
%     removes connected components (individual stream 'trees') that have
%     less or equal the number of channel heads 
%
%     'rmnodes' STREAMobj
%     removes nodes in S that belong another stream network S2.
%
%     'tributaryto' instance of STREAMobj
%     returns the stream network that is tributary to a stream (network) 
%
%     'tributaryto2' instance of STREAMobj
%     same as 'tributaryto' but tributaries include the pixel of the
%     receiving stream.
%
%     'shrinkfromtop' scalar distance
%     removes parts of the stream network that are within the specified
%     distance from the channelheads.
%
%     'interactive'  string
%        'polyselect': plots the stream network and starts a polygon tool to
%                      select the stream network of interest.
%        'reachselect': select a reach based on two locations on the network
%        'channelheadselect': select a number of channel heads and derive 
%                      stream network from them.
%
%
% Output arguments
%
%     S2    modified stream network (class: STREAMobj)
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S  = STREAMobj(FD,flowacc(FD)>1000);
%
%     % Extract the stream network that is tributary to the main trunk
%     St = modify(S,'tributaryto',trunk(S));
%     plot(S)
%     hold on
%     plot(St,'r')
%     holf off
%
%     % Split the stream network into the network downslope and upslope to
%     % the 1200 m contour
%     C = dilate(DEM-1200,ones(3)) > 0 & erode(DEM-1200,ones(3)) < 0;
%     C.Z = bwmorph(C.Z,'skel',inf);
%     C.Z = bwmorph(C.Z,'diag');
%     Su  = modify(S,'upstreamto',C);
%     Sl  = modify(S,'downstreamto',C);
%     plot(Su)
%     hold on
%     plot(Sl,'r')
%     hold off
%
%     
% See also: STREAMobj, STREAMobj/trunk, demo_modifystreamnet
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 16. January, 2017

narginchk(3,3)

% Parse Inputs
p = inputParser;         
p.FunctionName = 'modify';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));

addParamValue(p,'streamorder',[]);
addParamValue(p,'distance',[],@(x) isnumeric(x) && numel(x)<=2);
addParamValue(p,'maxdsdistance',[],@(x) isscalar(x) && x>0);
addParamValue(p,'interactive',[],@(x) ischar(validatestring(x,{'polyselect','reachselect','channelheadselect'})));
addParamValue(p,'tributaryto',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'tributaryto2',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'shrinkfromtop',[],@(x) isnumeric(x) && isscalar(x) && x>0);
addParamValue(p,'upstreamto',[],@(x) isa(x,'GRIDobj') || isnumeric(x));
addParamValue(p,'downstreamto',[],@(x) isa(x,'GRIDobj') || isnumeric(x));
addParamValue(p,'rmconncomps',[],@(x) isnumeric(x) && x>0 && isscalar(x));
addParamValue(p,'rmconncomps_ch',[],@(x) isnumeric(x) && x>=0 && isscalar(x));
addParamValue(p,'rmnodes',[],@(x) isa(x,'STREAMobj'));

parse(p,S,varargin{:});
S   = p.Results.S;


if ~isempty(p.Results.streamorder)
%% streamorder    
    so = streamorder(S);
    if isnumeric(p.Results.streamorder)
        validateattributes(p.Results.streamorder,{'numeric'},...
            {'scalar','integer','>',0,'<=',max(so)});
        I = so == p.Results.streamorder;
    else
        try
            % expecting relational operator
            sothres = nan;
            counter = 1;
            while isnan(sothres)
                sothres = str2double(p.Results.streamorder(counter+1:end));
                relop   = p.Results.streamorder(1:counter);
                counter = counter+1;
            end
        catch ME
            rethrow(ME)
        end

        relop = str2func(relop);
        I = relop(so,sothres);
        
    end


elseif ~isempty(p.Results.distance);
%% distance        
    d = S.distance;
    maxdist  = max(d);
    distrange = p.Results.distance;
    
    
    validateattributes(distrange,{'numeric'},...
        {'>=',0,'<=',maxdist});
    
    if isscalar(distrange);
        distrange(2) = maxdist;
    end
    
    if distrange(2) <= distrange(1)
        error('TopoToolbox:wronginput',...
            'The minimum distance must not be equal to or greater than the maximum distance');
    end
        
    
    I = d>=distrange(1) & d<=distrange(2);               % + norm([S.cellsize S.cellsize]);
    
elseif ~isempty(p.Results.maxdsdistance);
%% maximmum downstream distance    
    d = distance(S,'min_from_ch');
    I = d <= p.Results.maxdsdistance;
    
elseif ~isempty(p.Results.upstreamto);
%% upstream to    
    
    II = p.Results.upstreamto;
    
    if isa(II,'GRIDobj')
        validatealignment(S,II);
        II.Z = logical(II.Z);
    else
        % II contains linear indices into the DEM from which S was derived
        IX = II;
        I = ismember(IX,S.IXgrid);
        if ~all(I)
            error('TopoToolbox:modify','Some linear indices are not located on the stream network')            
        end
        II = GRIDobj(S,'logical');
        II.Z(IX) = true;
    end
  
    I = false(size(S.x));
    for r = numel(S.ix):-1:1;
        I(S.ix(r)) = II.Z(S.IXgrid(S.ixc(r))) || I(S.ixc(r));
    end
        
elseif ~isempty(p.Results.downstreamto);
%% downstream to    
    
    II = p.Results.downstreamto;
    
    if isa(II,'GRIDobj')
        validatealignment(S,II);
        II.Z = logical(II.Z);
    else
        % II contains linear indices into the DEM from which S was derived
        IX = II;
        I = ismember(IX,S.IXgrid);
        if ~all(I)
            error('TopoToolbox:modify','Some linear indices are not located on the stream network')            
        end
        II = GRIDobj(S,'logical');
        II.Z(IX) = true;
    end
    
    I = false(size(S.x));
    for r = 1:numel(S.ix);
        I(S.ixc(r)) = II.Z(S.IXgrid(S.ix(r))) || I(S.ix(r)) || I(S.ixc(r));
    end

elseif ~isempty(p.Results.tributaryto);
%% tributary to
    Strunk = p.Results.tributaryto;
    
    II = ismember(S.IXgrid,Strunk.IXgrid);
    I = false(size(S.x));
    for r = numel(S.ix):-1:1;
        I(S.ix(r)) = (II(S.ixc(r)) || I(S.ixc(r))) && ~II(S.ix(r));
    end

elseif ~isempty(p.Results.tributaryto2)
%% tributary to 2
    Strunk = p.Results.tributaryto2;
    [~,~,~,S] = intersectlocs(Strunk,S);
    return
    
elseif ~isempty(p.Results.shrinkfromtop)
%% shrink from top
    d = distance(S,'max_from_ch');
    I = d > p.Results.shrinkfromtop;
    
elseif ~isempty(p.Results.rmconncomps)
%% remove connected conn comps
    cc = conncomps(S);
    d  = S.distance;
    md = find(accumarray(cc,d,[max(cc) 1],@max) > p.Results.rmconncomps);
    I = ismember(cc,md);

elseif ~isempty(p.Results.rmconncomps_ch)
%% remove connected conn comps
    cc = conncomps(S);
    d  = distance(S,'nr_of_ch');
    md = find(accumarray(cc,d,[max(cc) 1],@max) > p.Results.rmconncomps_ch);
    I = ismember(cc,md);
    
elseif ~isempty(p.Results.rmnodes)
    I = ~ismember(S.IXgrid,p.Results.rmnodes.IXgrid);
    
elseif ~isempty(p.Results.interactive);
%% interactive    
    figure
    plot(S,'k'); axis image
    
    switch validatestring(p.Results.interactive,{'polyselect','reachselect','channelheadselect'})
        case 'channelheadselect'
            title('map channel heads and enter any key to finalize')
            hold on
            xy = streampoi(S,'channelhead','xy');
            scatter(xy(:,1),xy(:,2));
            hold off
            set(gcf,'WindowKeyPressFcn',@(k,l) uiresume);
            xys   = [];
            while true
                try
                    hpstart = impoint('PositionConstraintFcn',@getnearestchanhead);
                    
                    setColor(hpstart,[0 1 0])
%                     idstart = addNewPositionCallback(hpstart,@drawpath);
                    setPosition(hpstart,getPosition(hpstart))
                    xys = [xys; getPosition(hpstart)];
                catch
                    break
                end
            end
            
%             IX = coord2ind(S,xys(:,1),xys(:,2));
            I = ismember([S.x S.y],xys,'rows');
            
            for r = 1:numel(S.ix)
                I(S.ixc(r)) = I(S.ix(r)) || I(S.ixc(r));
            end
            
        case 'polyselect'
            title('create a polygon and double-click to finalize')
            hp = impoly;
            pos = wait(hp);
            pos(end+1,:) = pos(1,:);

            I = inpolygon(S.x,S.y,pos(:,1),pos(:,2));
            hold on
            plot(pos(:,1),pos(:,2),'-');
            hold off
            delete(hp);
        case 'reachselect'
            ix = [];
            hpath = [];
            ixpath = [];
            ixend  = 0;
            ixcix = zeros(numel(S.x),1);
            ixcix(S.ix) = 1:numel(S.ix);
            title('set upper (green) reach location')
            hpstart = impoint('PositionConstraintFcn',@getnearest);
            setColor(hpstart,[0 1 0])
            idstart = addNewPositionCallback(hpstart,@drawpath);
            setPosition(hpstart,getPosition(hpstart))
            
            title('set lower (red) reach location')
            hpend   = impoint('PositionConstraintFcn',@getnearestend);            
            setColor(hpend,  [1 0 0])
            idend = addNewPositionCallback(hpend,@drawpath);
            setPosition(hpend,getPosition(hpend))
            
            title('move reach locations and press any key to extract reach')
            set(gcf,'WindowKeyPressFcn',@(k,l) uiresume);
            uiwait
            
            I = false(size(S.x));
            I(ixpath) = true;
            
            delete(hpstart)
            delete(hpend);


    end
end

%% clean up
if exist('I','var');

L = I;
I = L(S.ixc) & L(S.ix);

S.ix  = S.ix(I);
S.ixc = S.ixc(I);

IX    = cumsum(L);
S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

S.x   = S.x(L);
S.y   = S.y(L);
S.IXgrid   = S.IXgrid(L);

% if ~isempty(p.Results.interactive) && p.Results.interactive;
%     hold on
%     plot(S,'r')
%     hold off
% end

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function posn = getnearest(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((S.x-pos(1)).^2 + (S.y-pos(2)).^2);
    posn= [S.x(ix) S.y(ix)];

end

function posn = getnearestchanhead(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((xy(:,1)-pos(1)).^2 + (xy(:,2)-pos(2)).^2);
    posn= xy(ix,:);

end

function posn = getnearestend(pos)
    % GETNEAREST  snap to nearest location on stream network
        [~,ixend] = min((S.x-pos(1)).^2 + (S.y-pos(2)).^2);
        posn= [S.x(ixend) S.y(ixend)];

end

function drawpath(pos)
    % DRAWPATH   
    
    % bring hpstart and hpend to top

    ixcc = ix;
    ixpath = ix;
    while ixcix(ixcc) ~= 0 && ixcc ~= ixend; %S.ixc(ixcix(ixcc)) ~= ixend;
        ixpath(end+1) = S.ixc(ixcix(ixcc));
        ixcc = ixpath(end);
    end

    if ishandle(hpath)
        set(hpath,'Xdata',S.x(ixpath),'Ydata',S.y(ixpath));
    else
        hold on
        hpath = plot(S.x(ixpath),S.y(ixpath),'r','LineWidth',1.5);
        hold off
    end
end 

        
end
            
            
        
    