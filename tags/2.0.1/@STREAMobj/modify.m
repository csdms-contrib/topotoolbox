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
%     different user-defined criteria.
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
%     let's you define a minimum and maximum distance from which the stream 
%     network starts and ends based on the distances given in S. A scalar is
%     interpreted as minimum upstream distance and a two element vector is
%     interpreted as minimum and maximum upstream distance to which the
%     stream network is cut.
%
%     'upstreamto' logical GRIDobj
%     returns the stream network upstream to true pixels in the logical
%     raster of GRIDobj. Note that, if the grid contains linear features
%     (e.g. a fault), the line should be 4 connected. Use
%     bwmorph(I.Z,'diag') to establish 4-connectivity in a logical raster
%     I.
%
%     'downstreamto' logical GRIDobj 
%     returns the stream network downstream to true pixels in the logical
%     raster of GRIDobj. Note that, if the grid contains linear features
%     (e.g. a fault), the foreground (true pixels) should be 4 connected.
%     Use bwmorph(I.Z,'diag') to establish 4-connectivity in a logical
%     raster I.
%
%     'tributaryto' instance of STREAMobj
%     returns the stream network that is tributary to a stream (network) 
%
%     'interactive' string
%     'polyselect': plots the stream network and starts a polygon tool to
%     select the stream network of interest.
%     'reachselect': select a reach based on two locations on the network
%
%
% Output arguments
%
%     S2    modified stream network (class: STREAMobj)
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'mex',true,'preprocess','carve');
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
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. October, 2013

narginchk(3,3)

% Parse Inputs
p = inputParser;         
p.FunctionName = 'modify';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));

addParamValue(p,'streamorder',[]);
addParamValue(p,'distance',[],@(x) isnumeric(x) && numel(x)<=2);
addParamValue(p,'maxdistance',[],@(x) isscalar(x));
addParamValue(p,'interactive',[],@(x) ischar(validatestring(x,{'polyselect','reachselect'})));
addParamValue(p,'tributaryto',[],@(x) isa(x,'STREAMobj'));
addParamValue(p,'upstreamto',[],@(x) isa(x,'GRIDobj'));
addParamValue(p,'downstreamto',[],@(x) isa(x,'GRIDobj'));

parse(p,S,varargin{:});
S   = p.Results.S;

   
if ~isempty(p.Results.streamorder)
    
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
            while isnan(sothres);
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
    
    distance = S.distance;
    maxdist  = max(distance);
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
        
    
    I = distance>=distrange(1) & distance<=distrange(2);               % + norm([S.cellsize S.cellsize]);
    
elseif ~isempty(p.Results.upstreamto);
    
    II = p.Results.upstreamto;
    validateattributes(II,{'GRIDobj'},{'scalar'});
    validatealignment(S,II);

    II.Z = logical(II.Z);
    I = false(size(S.x));
    for r = numel(S.ix):-1:1;
        I(S.ix(r)) = II.Z(S.IXgrid(S.ixc(r))) || I(S.ixc(r));
    end
        
elseif ~isempty(p.Results.downstreamto);
    
    II = p.Results.downstreamto;
    validateattributes(II,{'GRIDobj'},{'scalar'});
    validatealignment(S,II);

    II.Z = logical(II.Z);
    I = false(size(S.x));
    for r = 1:numel(S.ix);
        I(S.ixc(r)) = II.Z(S.IXgrid(S.ix(r))) || I(S.ix(r)) || I(S.ixc(r));
    end

elseif ~isempty(p.Results.tributaryto);
    Strunk = p.Results.tributaryto;
    
    II = ismember(S.IXgrid,Strunk.IXgrid);
    I = false(size(S.x));
    for r = numel(S.ix):-1:1;
        I(S.ix(r)) = (II(S.ixc(r)) || I(S.ixc(r))) && ~II(S.ix(r));
    end
    
elseif ~isempty(p.Results.interactive);
    
    figure
    plot(S,'k'); axis image
    
    switch validatestring(p.Results.interactive,{'polyselect','reachselect'})
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

function posn = getnearest(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((S.x-pos(1)).^2 + (S.y-pos(2)).^2);
    posn= [S.x(ix) S.y(ix)];

end

function posn = getnearestend(pos)
    % GETNEAREST  snap to nearest location on stream network
        [~,ixend] = min((S.x-pos(1)).^2 + (S.y-pos(2)).^2);
        posn= [S.x(ixend) S.y(ixend)];

end

function drawpath(pos)
    % DRAWPATH   
    ixcc = ix;
    ixpath = [];
    while ixcix(ixcc) ~= 0 && S.ixc(ixcix(ixcc)) ~= ixend;
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
            
            
        
    