function tf = isequal(S1,S2,type)

%ISEQUAL Determine whether two STREAMobjs are equal
%
% Syntax 
%
%     tf = isequal(S1,S2)
%     tf = isequal(S1,S2,type)
%
% Description
%
%     isequal determines whether two instances of the class STREAMobj are
%     equal. Equality means that both have the same referencing, the same
%     geometry, and topology. 
%
% Input arguments
%
%     S1, S2   instances of STREAMobj
%     type     'geometry' or {'topology'}. 'topology' requires that both
%              STREAMobjs have the same network structure, including the
%              same sorting. 'geometry' requires only that that the two 
%              networks have the same nodes.
%      
% Output arguments
%
%     tf       true or false (scalar)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S2 = trunk(S);
%     isequal(S,S2)    
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. January, 2020

if nargin == 2
    type = 'topology';
else    
    type = validatestring(type,{'referencing','geometry','topology'});
end


% georeferencing and refmat
tf = isequal(S1.refmat,S2.refmat);
if ~tf || strcmp(type,'referencing'); return; end

% geometry
tf = isequal(S1.IXgrid,S2.IXgrid);
if ~tf || strcmp(type,'geometry'); return; end

% topology
tf = isequal(S1.ix,S2.ix) && isequal(S1.ixc,S2.ixc);
if ~tf; return; end