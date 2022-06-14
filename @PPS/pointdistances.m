function d = pointdistances(P,varargin)

%POINTDISTANCE inter-point distance
%
% Syntax
%
%     d = pointdistances(P)
%     d = pointdistances(P,pn,pv,...)
%
% Description
%
%     pointdistances calculates inter-point distances on stream networks
%     using different methods.
% 
% Input arguments
%
%     P       instance of PPS
%     
%     Parameter name/value pairs
%     
%     'type'        {'graph'},'digraph','nearest'
%     'direction'   'up', 'down' (only applicable for 'digraph' and 
%                   'nearest')
%                   'both' (only applicable for 'nearest')
%     'd3d'         false, if distance in 2d (default) or true, if in 3d 
%     'val'         By default, pointdistances measures the distance along
%                   the stream network. However, you can also other
%                   node-attribute lists (e.g. elevation)
%     'output'      {'matrix'} or 'vector'. If 'type', 'nearest' then the 
%                   distance is always returned as vector. Note that if you
%                   use 'type' = 'digraph', then the matrix is asymmetric.
%                   Don't use 'vector' output in this case. 
%
% Output arguments
%
%     d       matrix or vector of distances
%
% Example
%
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = removeshortstreams(S,100);
%     S = clean(S);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     d = pointdistances(P);
%     imagesc(d);
%     axis image
%     title('Distance matrix');
%
%
% See also: PPS, PPS/netdist
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. June, 2020


% Check input arguments
p = inputParser;
p.FunctionName = 'PPS/pointdistances';
addParameter(p,'output','matrix');
addParameter(p,'type','graph');
addParameter(p,'d3d',false);
addParameter(p,'extendednetwork',false);
addParameter(p,'val',[]);
addParameter(p,'direction','down');
addParameter(p,'sourcepoints',[]);
addParameter(p,'targetpoints',[]);
% Parse
parse(p,varargin{:});

ix  = P.S.ix;
ixc = P.S.ixc;
d   = nodedistance(P,'d3d',p.Results.d3d,'val',p.Results.val);

switch lower(p.Results.type)
    case {'graph','digraph'}
        switch lower(p.Results.type)
            case 'graph'
                G = graph(ix,ixc,d);
            case 'digraph'
                switch p.Results.direction
                    case 'down'
                        G = digraph(ix,ixc,d);
                    case 'up'
                        G = digraph(ixc,ix,d);
                end
        end
        
        tf = hasduplicates(P);
        if tf            
            % sigh... has duplicates 
            if ~p.Results.extendednetwork
                [P2,~,locb] = unique(P.PP,'stable');
                d2 = distances(G,P2,P2,'Method','positive');
                d  = d2(locb,locb);
            else           
                [G,pp,pback] = extendednetwork(P);
                d = distances(G,pp,pp,'Method','positive');
                d(pback,pback) = d;
            end
         
        else
            % no duplicates
            d = distances(G,P.PP,P.PP,'Method','positive');
        end
        
    case 'nearest'
        
        switch p.Results.direction
            case 'down'
                
                dm = nan(size(P.S.IXgrid));
                pp = ~points(P,'nal');
                for r = numel(ix):-1:1
                    if pp(ixc(r))
                        dm(ix(r)) = dm(ixc(r)) + d(r);
                    else
                        dm(ix(r)) = d(r);
                    end
                end
                d  = dm(P.PP);
                return
            case 'up'
                
                dm = nan(size(P.S.IXgrid));
                pp = ~points(P,'nal');
                for r = 1:numel(ix)
                    if pp(ix(r))
                        dm(ixc(r)) = min(dm(ix(r)) + d(r),dm(ixc(r)));
                    else
                        dm(ixc(r)) = d(r);
                    end
                end
                d  = dm(P.PP);
                return
                
            case 'both'
                [P.PP,~,locb] = unique(P.PP,'stable');
                [v,d] = voronoi(P);
                % get points at
                II = v(P.S.ix) ~= v(P.S.ixc);
                I  = false(size(P.S.IXgrid));
                I(P.S.ix(II))  = true;
                I(P.S.ixc(II)) = true;
                
                d(~I) = inf;
                dd    = d(P.S.ix(II)) - d(P.S.ixc(II)); 
                
                
                
                d(P.S.ix(II))  = d(P.S.ix(II)) + dd/2;  
                d(P.S.ixc(II)) = d(P.S.ixc(II)) - dd/2;  
                
                
                d  = accumarray(v,d,[npoints(P) 1],@min);
                d  = d*2;
                d  = d(locb);
                
                return
                
        end
end

outp = validatestring(p.Results.output,{'matrix', 'vector', 'edgelist'});

switch outp
    case 'matrix'
    case 'vector'
        d = squareform(d,'tovector')';
    case 'edgelist'
        [i,j,d] = find(tril(d));
        I = isfinite(d);
        d = [P.PP(i(I)) P.PP(j(I)) d(I)];
end
