function d = distance(S,type)

%DISTANCE return node attribute list with distances along the stream network
%
% Syntax
%
%     d = distance(S,type)
%     d = distance(S,S2)
%
% Description
%
%     This function returns a node attribute list with a distance value for  
%     each node in the stream network S. Since distance can be calculated
%     in different ways for a channelnetwork, the second input argument
%     allows choosing between a number of different distance metrics. If
%     the second input argument S2 is a STREAMobj, distances will be based
%     on the distances of S2. Note that S must be a subset of S2 then! 
%
% Input arguments
%
%     S       STREAMobj
%     type    'from_outlet'      distance in upstream direction
%             'min_from_ch'      shortest from channelhead
%             'max_from_ch'      longest from channelhead
%             'mean_from_ch'     mean distance from channelheads
%             'nr_of_ch'         number of channelheads (not a distance
%                                measure)
%             'node_to_node'     distance between each node and its
%                                downstream neighbor
%             'accumdownstream'  downstream accumulated distance
%             'from_trunk'       distance in upstream direction from trunk
%                                stream
%     S2      STREAMobj of which S is a subset.
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. September, 2016

if nargin == 1;
    d = S.distance;
    return
end

if ischar(type)
    validtypes = {...
        'from_outlet', ... % distance in upstream direction
        'min_from_ch',... % shortest from channelhead
        'max_from_ch',... % longest from channelhead
        'mean_from_ch',... % mean distance from channelheads
        'nr_of_ch'...
        'node_to_node'...
        'accumdownstream' ...
        'from_trunk' ...
        };
    
    type = validatestring(type,validtypes,'distance','type',2);
    
    % handle the simple case which is implement as dynamic property of S
    switch type
        case 'from_outlet'
            d = S.distance;
            return;
    end
    
    % preallocate array
    d = nan(size(S.x));
    
    % distance between two nodes
    dedge = sqrt((S.x(S.ixc)-S.x(S.ix)).^2 + (S.y(S.ixc)-S.y(S.ix)).^2);
    
    switch type
        case 'node_to_node'
            d(:) = 0;
            d(S.ix) = dedge;
            return
    end
    
    % identify channel heads
    M  = sparse(double(S.ix),double(S.ixc),true,numel(d),numel(d));
    CH = (sum(M,1)'==0) & (sum(M,2) ~= 0);
    
    % set distance at channel heads to zero
    d(CH) = 0;
    switch type
        case 'min_from_ch'
            for r = 1:numel(S.ix);
                d(S.ixc(r)) = min(d(S.ix(r)) + dedge(r),d(S.ixc(r)));
            end
        case 'max_from_ch'
            for r = 1:numel(S.ix);
                d(S.ixc(r)) = max(d(S.ix(r)) + dedge(r),d(S.ixc(r)));
            end
        case 'nr_of_ch'
            nrc = numel(d);
            d = full((speye(nrc)-M')\double(CH));
        case 'mean_from_ch'
            nrc = numel(d);
            nr = (speye(nrc)-M')\double(CH);
            
            d  = accumarray(S.ixc,dedge,[nrc,1],@mean,0).*nr;
            
            d  = (speye(nrc)-M')\d;
            d  = full(d./nr);
        case 'accumdownstream'
            d(:) = 0;
            for r = 1:numel(S.ix);
                d(S.ixc(r)) = d(S.ix(r))+dedge(r) + d(S.ixc(r));
            end
        case 'from_trunk'
            d(:) = 0;
            St   = trunk(S);
            I    = ~ismember(S.IXgrid,St.IXgrid);
            for r = numel(S.ix):-1:1
                if I(S.ix(r))
                    d(S.ix(r)) = d(S.ixc(r)) + dedge(r);
                end
            end
                    
    end
    
else
    S2 = type;
    d2 = S2.distance;
    [I,inS2] = ismember(S.IXgrid,S2.IXgrid);
    if ~all(I)
        error('S is not a subset of S2');
    end
    d  = d2(inS2);
end
    
