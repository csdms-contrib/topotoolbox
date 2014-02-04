function OUT = streamorder(FD,WW,type)


% calculate Strahler Stream Order Grid from FLOWobj
%
% Syntax
%
%     S = streamorder(FD,W)
%     S = streamorder(FD,W,type)
%
% Description
%
%     Stream order is a way to classify rivers based on a hierarchy of
%     tributaries. In the Strahler ordering scheme, first-order streams
%     don't have tributaries. When two first-order streams come together,
%     they form a second-order stream. When two second-order streams come
%     together, they form a third-order stream. When two streams of
%     different order confluence, they form a stream of maximum order of
%     both.
%
%     In the Shreve ordering scheme,a ll links with no tributaries are
%     assigned a magnitude (order) of one. Magnitudes are additive
%     downslope. When two links intersect, their magnitudes are added and
%     assigned to the downslope link.
%
%     streamorder returns the Strahler or Shreve Stream Order based on an
%     instance of the FLOWobj (FD) and a stream grid (W). W is either a
%     logical matrix or a GRIDobj that contains a logical matrix, that has
%     ones where streams are and zeros elsewhere.
% 
%     The output S is a GRIDobj and contains the Strahler or Shreve Order
%     for each cell. Non-channel cells are set to zero.
%
%     Stream order can also be calculated based on an instance of
%     STREAMobj.
%
% Input
%   
%     FD        FLOWobj
%     W         channel grid (logical matrix, true where channels/
%               channelheads are, or GRIDobj)
%     type      'strahler' (default) or 'shreve'
%
% Output
%
%     S         stream order grid (GRIDobj)
%
% Example:
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S = streamorder(FD,flowacc(FD)>1000);
%     imagesc(S)
% 
% See also: FLOWobj, FLOWobj/flowacc, STREAMobj/streamorder
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013


%% check input arguments
narginchk(2,3)
validatealignment(FD,WW);

if nargin == 2;
    type = 'strahler';
else
    type = validatestring(type,{'strahler','shreve'},'FLOWobj/streamorder','type',3);
end

if ~strcmpi(FD.type,'single')
    error('TopoToolbox:FLOWobj','streamorder requires FLOWobj.type to be single')
end

if isa(WW,'GRIDobj');
    W   = WW.Z;
else
    W   = WW;
end

switch lower(type)
    case 'strahler'
        % Strahler stream order
        
        S   = uint16(W);
        offset = uint16(1);
        VIS    = false(FD.size);
        
        for r = 1:numel(FD.ix);
            if W(FD.ix(r))
                if (S(FD.ixc(r)) == S(FD.ix(r))) && VIS(FD.ixc(r));
                    S(FD.ixc(r)) = S(FD.ixc(r))+offset;
                    VIS(FD.ixc(r)) = true;
                else
                    S(FD.ixc(r)) = max(S(FD.ix(r)),S(FD.ixc(r)));
                    VIS(FD.ixc(r)) = true;
                end
            end
        end
        
    case 'shreve'
        S = streampoi(FD,W,'channelheads');
        S = flowacc(FD,S);
        S = uint16(S.Z);        
end



%% Prepare Output
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = S;
OUT.zunit = '';
OUT.name  = [type ' stream order'];
