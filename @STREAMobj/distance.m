function d = distance(S,type)



validtypes = {...
              'from_outlet', ... % distance in upstream direction
              'min_from_ch',... % shortest from channelhead
              'max_from_ch',... % longest from channelhead
              'mean_from_ch',... % mean distance from channelheads
              'nr_of_ch'...
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
end
        
        
        
        
        
        