function val = info(S,m)

%INFO meta information about STREAMobj
%
% Syntax
%
%     in = info(S)
%
% Description
%
%     info(S) summarizes meta information about a stream network in a 
%     structure array.
%
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016


validargs = {'nrnodes' ...
             'boundingbox' ...
             'centroid', ...
             'totallength' ...
             'nredges' ...
             'nrchannelheads' ...
             'nroutlets' ...
             'nrconfluences' ...
             'maxstreamorder' ...
             'maxdistance'};
if nargin == 1
    for r = 1: numel(validargs)
        val.(validargs{r}) = info(S,validargs{r});
    end
    return
else
    m = validatestring(m,validargs,'STREAMobj/info','type',2);
end
    

switch m
    case 'nrnodes'
        val = numel(S.IXgrid);
    case 'boundingbox'
        val = [min(S.x) max(S.x) min(S.y) max(S.y)];
    case 'centroid'
        val = mean([S.x S.y]);
    case 'totallength'
        d   = distance(S,'node_to_node');%sqrt((S.x(S.ix)-S.x(S.ixc)).^2 + (S.y(S.ix)-S.y(S.ixc)).^2);
        val = sum(d);
    case 'nredges'
        val = numel(S.ix);
    case 'nrchannelheads'
        nrc = info(S,'nrnodes');
        M = sparse(S.ix,S.ixc,true,nrc,nrc);
        val = full(sum(sum(M,1) == 0));
    case 'nroutlets'
        nrc = info(S,'nrnodes');
        M = sparse(S.ix,S.ixc,true,nrc,nrc);
        val = full(sum(sum(M,2) == 0));
    case 'nrconfluences'
        nrc = info(S,'nrnodes');
        M = sparse(S.ix,S.ixc,true,nrc,nrc);
        val = full(sum(sum(M,1) > 1));
    case 'maxstreamorder'
        s   = streamorder(S,'strahler');
        val = max(s);
    case 'maxdistance'
        val = max(S.distance);
end
