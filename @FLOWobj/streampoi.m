function OUT = streampoi(FD,WW,type,outformat)

% stream points of interest
%
% Syntax
%
%     V = streampoi(FD,W,type)
%     V = streampoi(FD,W,type,outformat)
%
% Description
%
%     streampoi returns points of interest in a stream network. These
%     points can be 'channelheads', 'confluences' or 'outlets'.
%
% Input arguments
%   
%     FD          flow direction (FLOWobj)
%     W           channel grid (logical matrix, true where channels/
%                 channelheads are, or GRIDobj) 
%     type        'channelheads' (default), 'confluences' or 'outlets'
%     outformat   'GRIDobj' (default): returns an instance of GRIDobj 
%                 'xy': nx2 coordinate matrix with x and y coordinates of 
%                       n points 
%                 'rc': nx2 matrix with row and column subscripts
%                 'ix': nx1 vector with linear indices
%
% Output
%
%     V         output as specified in the option 'outformat'
%
% Example
% 
% See also: FLOWobj, FLOWobj/flowacc, FLOWobj, STREAMobj/streampoi
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013


%% check input arguments
narginchk(2,4)
validatealignment(FD,WW);

if nargin == 2;
    type = 'channelheads';
    outformat = 'GRIDobj';
else 
    type = validatestring(type,{'channelheads','confluences','outlets'},'streampoi','type',3);
    if nargin > 3
        outformat = validatestring(outformat,{'gridobj','xy','rc','ix'},'streampoi','outformat',4);
    else       
        outformat = 'GRIDobj';
    end  
end

if ~strcmpi(FD.type,'single')
    error('TopoToolbox:FLOWobj','streamorder requires FLOWobj.type to be single')
end

if isa(WW,'GRIDobj');
    W   = WW.Z;
else
    W   = WW;
end

I   = W(FD.ix);
ix  = FD.ix(I);
ixc = FD.ixc(I);

siz = FD.size;
nrc = prod(siz);

M  = sparse(double(ix),double(ixc),true,nrc,nrc);

switch type
    case 'channelheads'
        
        V = (sum(M,1)'==0)& W(:);
        V = reshape(full(V),siz);
        
    case 'confluences'
        
        V = reshape(full(sum(M,1)'>1),siz);

    case 'outlets'

        V = reshape(full(sum(M,2)==0),siz) & W;
end


switch lower(outformat)
    case 'gridobj'
        %% Prepare Output
        % write output to GRIDobj
        OUT = copy2GRIDobj(FD);
        OUT.Z = V;
        OUT.zunit = '';
        OUT.name  = type;
    case 'xy'
        TEMP = copy2GRIDobj(FD);
        [OUT(:,1) OUT(:,2)] = ind2coord(TEMP,find(V));
    case 'rc'
        [OUT(:,1) OUT(:,2)] = find(V);
    case 'ix'
        OUT = find(V);
end
        
        
        

