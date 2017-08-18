function varargout = intersectlocs(S1,S2)

%INTERSECTLOCS Derive locations where two STREAMobj start to have a common network
%
% Syntax
%
%     IX = intersectlocs(S,Sl)
%     [IX,x,y,Slc] = intersectlocs(S,Sl)
%
% Description
%
%     intersectlocs returns for two stream networks S and S1 (STREAMobj)  
%     derived from a FLOWobj the locations where S1 starts to have common
%     flowpath with S. The function may be used to map confluences of a low
%     order stream network with a  network that contains only higher order
%     streams. It may also be useful to map intersections of landslide
%     runout zones with a river network.
%
%     Note that intersectlocs might return unexpected results if both 
%     STREAMobjs were not derived from the same FLOWobj.
%
% Input arguments
%
%     S     STREAMobj
%     S1    STREAMobj
%     
% Output arguments
%
%     IX    Linear indices of locations. Indices refer the DEM from 
%           which FLOWobj and STREAMobj where derived.
%           where derived.
%     x,y   coordinate pairs of intersections
%     Slc   stream network same as S1 but comprises only the streams
%           between channel heads and intersection locations.
%
% Example
%
%     % Map locations where landslides hit rivers
%     % http://topotoolbox.wordpress.com/2014/09/09/landslides-hit-rivers/
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c','mex',true);
%     S = STREAMobj(FD,flowacc(FD)>1000);
%     LS = GRIDobj(DEM);
%     LS.Z = logical(LS.Z);
%     LS.Z(randperm(numel(LS.Z),100)) = true;
%     LSi = influencemap(FD,LS);
%     Sls = STREAMobj(FD,LSi);
%     [IX,x,y,Slc] = intersectlocs(S,Sls);
%     subplot(2,1,1)
%     imageschs(DEM);
%     hold on
%     plot(S,'k');
%     plot(x,y,'*w');
%     plot(Slc,'w');
% 
%     subplot(2,1,2)
%     plotdz(S,DEM,'annotation',IX);
%
%
% See also: STREAMobj, STREAMobj/intersect
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. September, 2014



% identify the nodes that both networks have in common
SharedNetworkNodes = ismember(S2.IXgrid,S1.IXgrid);

% identify edges in S2 whose whose receiver nodes do belong to the 
% common nodes and whose giver nodes do not.
I  = SharedNetworkNodes(S2.ixc) & ~SharedNetworkNodes(S2.ix);
 
varargout{1} = S2.IXgrid(S2.ixc(I));
if nargout > 1
    varargout{2} = S2.x(S2.ixc(I));
    varargout{3} = S2.y(S2.ixc(I));
    
    II = false(size(S2.x));
    II(S2.ixc(I)) = true;
    % get network upstream to locations
    for r = numel(S2.ix):-1:1;
        if II(S2.ixc(r))
            II(S2.ix(r)) = II(S2.ixc(r)) & ~SharedNetworkNodes(S2.ix(r));
        end
    end
    if ~isempty(II);

        L = II;
        I = L(S2.ixc) & L(S2.ix);

        S2.ix  = S2.ix(I);
        S2.ixc = S2.ixc(I);
    
        IX    = cumsum(L);
        S2.ix  = IX(S2.ix);
        S2.ixc = IX(S2.ixc);

        S2.x   = S2.x(L);
        S2.y   = S2.y(L);
        S2.IXgrid   = S2.IXgrid(L);
    
    end
    varargout{4} = S2;
end