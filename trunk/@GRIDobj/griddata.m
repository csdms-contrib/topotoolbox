function Z = griddata(x,y,z,DEM,varargin)

% Use different techniques to spatially interpolate values at grid locations
%
% Syntax
%
%     Z = griddata(x,y,z,DEM,method)
%
% Description
%
%     griddata overloads various methods that come with core Matlab and
%     adds a few other methods to obtain spatially predictions at grid
%     locations based on scattered data. Currently supported methods are
%
%
% 
% use scatteredInterpolant in future versions



Z = DEM;
Z.name = ['griddata (' p.Results.method ')'];

switch lower(method)
    case 'aggregate'
        % requires 'aggfun' as additional input
        fun = p.Results.aggfun;
        if isa(fun,'function_handle');
        else
            fun = str2func(fun);
        end
        
        if isinteger(z)
            nanval = zeros(1,1,class(z));
        elseif islogical(z)
            nanval = false;
        else
            nanval = nan;
        end
        
        IX = coord2ind(DEM,x,y);
        Z.Z = reshape(accumarray(IX,z,[prod(DEM.size) 1],fun),DEM.size,nanval);
    case 'idw'
        % requires 'beta','searchradius' and 'nrneighs' as additional
        % inputs
        [X,Y] = getcoordinates(DEM);
        Z.Z   = blockproc(DEM.Z,bestblk(DEM.size,500),@idw,'UseParallel',true);
end

        
        
        



    function OUT = idw(B)

% block_struct.blockSize: A two-element vector, [rows cols], that specifies
% the size of the block data. If a border has been specified, the size does
% not include the border pixels.
%
% block_struct.data: M-by-N or M-by-N-by-P matrix of block data
%
% block_struct.imageSize: A two-element vector, [rows cols], that specifies
% the full size of the input image.
%
% block_struct.location: A two-element vector, [row col]
        Xi = X(B.location(2):B.location(2)+B.blockSize(2)-1);
        Yi = Y(B.location(1):B.location(1)+B.blockSize(1)-1);
        
        [Xi,Yi] = meshgrid(Xi,Yi);
        
        [IX,D] = knnsearch([x(:) y(:)],[Xi(:) Yi(:)],'K',p.Results.nrneigh);
        D(D>p.Results.searchradius) = inf;

        GAMMA = 1./(D.^p.Results.beta);
        sGAMMA = sum(GAMMA,2);
        GAMMA = bsxfun(@rdivide,GAMMA,sGAMMA);
        GAMMA(sGAMMA == 0) = nan;


        OUT = sum(z(IX).*GAMMA,2);
        OUT = reshape(OUT,B.blockSize);
        
        
        
    end

end
        
        


