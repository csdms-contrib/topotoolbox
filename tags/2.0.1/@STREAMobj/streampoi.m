function OUT = streampoi(S,type,outformat)

% stream points of interest
%
% Syntax
%
%     V = streampoi(S,type)
%     V = streampoi(S,type,outformat)
%
% Description
%
%     streampoi returns points of interest in a stream network. These
%     points can be 'channelheads', 'confluences' or 'outlets'.
%
% Input arguments
%   
%     FD          stream network (class: STREAMobj)
%     type        'channelheads' (default), 'confluences' or 'outlets'
%     outformat   'xy': nx2 coordinate matrix with x and y coordinates of 
%                       n points (default)
%                 'ix': nx1 vector with linear indices into an instance of
%                       GRIDobj with the same dimension as the GRIDobj from
%                       which S was derived.
%                 'logical': node attribute list (logical)
%
% Output
%
%     V           output as specified in outformat
%
% Example
% 
% See also: STREAMobj, FLOWobj/streampoi
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013


%% check input arguments
narginchk(1,3)

if nargin == 1;
    type = 'channelheads';
    outformat = 'xy';
else 
    type = validatestring(type,{'channelheads','confluences','outlets'},'streampoi','type',2);
    if nargin > 2
        outformat = validatestring(outformat,{'xy','ix','logical'},'streampoi','outformat',3);
    else       
        outformat = 'xy';
    end  
end

ix  = S.ix;
ixc = S.ixc;

nrc = numel(S.x);

M  = sparse(double(ix),double(ixc),true,nrc,nrc);

switch type
    case 'channelheads'
        
        V = (sum(M,1)'==0) & (sum(M,2) ~= 0);
        
    case 'confluences'
        
        V = sum(M,1)'>1;

    case 'outlets'

        V = (sum(M,2)==0) & (sum(M,1)'~=0);
end


switch outformat
    case 'xy'      
        ix   = find(V);
        ix   = ix(:);
        OUT  = [S.x(ix) S.y(ix)];
        
    case 'ix'
        ix   = find(V);
        ix   = ix(:);
        OUT  = S.IXgrid(ix);
        
    case 'logical'
        
        OUT  = full(V(:));
end