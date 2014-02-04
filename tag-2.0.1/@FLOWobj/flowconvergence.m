function N = flowconvergence(FD)

% compute flow convergence of a digital elevation model
% 
% Syntax
%
%     N = flowconvergence(FD)
%
% Description
%
%     flowconvergence computes the number of neighboring cells that 
%     drain into each cell. Values in N range between 0 on ridges to 8 in
%     pits. N can thus be used as a measure of flow convergence in a
%     digital elevation model. 
%
% Input
%
%     FD      flow direction (class: FLOWobj)
%
% Output
%
%     N       grid containing the sum of contributing neighbor cells. 
%             (class: GRIDobj)
%
% Example
% 
% 
%
% See also: FLOWOBJ
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. February, 2013


nr = zeros(FD.size,'uint8');
for r = 1:numel(FD.ix);
    nr(FD.ixc(r)) = nr(FD.ixc(r))+1;
end

% Prepare output
N = copy2GRIDobj(FD);
N.Z = nr;
N.name = 'flow convergence';