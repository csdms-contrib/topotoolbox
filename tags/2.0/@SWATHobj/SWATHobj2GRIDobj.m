function [OUT] = SWATHobj2GRIDobj(SW,DEM,type)
% convert SWATHobj to GRIDobj with swath-specific information
%
% Syntax
%
%     OUT = SWATHobj2GRIDobj(SW,DEM,type)
%
% Description
%
%     SWATHobj2GRIDobj(SW,DEM,type) 
%
% Input arguments
%
%     SW     instance of SWATHobj
%     DEM    instance of GRIDobj
%     type   string that specifies the type of output. Valid choices are
%            'ix', 'distx', 'disty', and 'z'. 'ix' writes down the index of
%            each item in the SWATHobj. 'distx' and 'disty' give the along
%            and across swath distance, respectively. 'z' writes down the
%            z-value of the accompanying GRIDobj.
%
% Output arguments
%
%     OUT    instance of GRIDobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>1000);
%     S = klargestconncomps(S,1);
%     S = removeshortstreams(S,1e3);
%     SW = SWATHobj(DEM,S,'smooth',300,'plot',false);
%     SW = tidyswath(SW,FD,'both');
%     DX = SWATHobj2GRIDobj(SW,DEM,'distx');
%     figure,imagesc(DX), hold on
%     plot(SW,'points',false,'labeldist',1e3), hold off
%     colorbar
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: November, 2013



narginchk(2,3)

if nargin == 2;
    type = 'z';
else
    type = validatestring(type,{'z','distx','disty','ix'});
end

% output grid
OUT = DEM;
OUT.Z(:) = nan;

for i = 1 : length(SW.Z)
    
    % delineate swath profile
    IM = true(size(SW.Z{i}));
    IM(isnan(SW.Z{i})) = false;
    B = bwboundaries(IM,4);
    A = false(DEM.size);
    for k = 1 : length(B)
        ix_bdy = sub2ind(size(IM),B{k}(:,1),B{k}(:,2));
        x_bdy = SW.X{i}(ix_bdy);
        y_bdy = SW.Y{i}(ix_bdy);
        if length(x_bdy)>2
            [newx,newy] = eqdistline(x_bdy,y_bdy,DEM.cellsize/3);
        else
            newx = x_bdy; newy = y_bdy;
        end
        ix_edge = coord2ind(DEM,newx,newy);
        ix_edge = unique(ix_edge);
        A(ix_edge) = true;
    end
    A = imfill(A,4,'holes');
    
    % swath locations
    ix_sw = find(~isnan(SW.Z{i}));
    x = SW.X{i}(ix_sw);
    y = SW.Y{i}(ix_sw);
    
    % grid locations
    ix_grid = find(A==1);
    [~,xmap,ymap] = GRIDobj2mat(DEM);
    [X,Y] = meshgrid(xmap,ymap);
    xq = X(ix_grid);
    yq = Y(ix_grid);
    [xq_edges,yq_edges] = ind2coord(DEM,ix_edge);
    
    switch type
        
        case 'z'
            OUT.Z(ix_grid) = DEM.Z(ix_grid);
            
        case 'distx' % longitudinal distance
            Dx = repmat(SW.distx{i}',length(SW.disty{i}),1);
            v = Dx(ix_sw);
            F = TriScatteredInterp(x,y,v);
            vq = F(xq,yq);
            OUT.Z(ix_grid) = vq;
            F = TriScatteredInterp(x,y,v,'nearest');
            vq = F(xq_edges,yq_edges);
            OUT.Z(ix_edge) = vq;
            
        case 'disty' % transverse distance
            Dy = repmat(SW.disty{i},1,length(SW.distx{i}));
            v = Dy(ix_sw);
            F = TriScatteredInterp(x,y,v);
            vq = F(xq,yq);
            OUT.Z(ix_grid) = vq;
            F = TriScatteredInterp(x,y,v,'nearest');
            vq = F(xq_edges,yq_edges);
            OUT.Z(ix_edge) = vq;
            
        case 'ix'
            OUT.Z(ix_grid) = i;
    
    end
    
end


    
    
    
    
    
    
    
    
    
    