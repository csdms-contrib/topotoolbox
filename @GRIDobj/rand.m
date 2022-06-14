function Z = rand(DEM,dist,params)

%RAND Compute a GRIDobj with random numbers
%
% Syntax
%
%      Z = rand(A)
%      Z = rand(A,'unif',[ll ul])
%      Z = rand(A,'normal', [mu sigma])
%      Z = rand(A,'randperm')
%      Z = rand(A,'randi',[ll ul])
%      Z = rand(A,pd)
%
% Description
%
%      rand returns a GRIDobj of random numbers. The raster has the same
%      dimensions and coordinates as the GRIDobj A. 
%    
%      rand(A,'unif',[ll ul]) returns uniform distributed numbers between
%      ll and ul.
%
%      rand(A,'normal', [mu sigma]) returns normal distributed random 
%      numbers with average mu and standard deviation sigma.
%
%      rand(A,'randperm') return a random permutation of integers in the 
%      range of 1:prod(A.size).
%
%      rand(A,'randi',[ll ul]) returns random integer values between ll and
%      ul.
%     
%      rand(A,pd) returns random numbers from the distribution defined by
%      pd (see makedist and fitdist).
%
% Input arguments
%
%      A           GRIDobj
%      ll,ul       lower and upper limits of random numbers.
%
% Output arguments
%
%      Z           GRIDobj
%
%
% See also: GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 31. May, 2022



if nargin == 1
    Z = GRIDobj(DEM);
    Z.Z = rand(DEM.size);
else
    if ischar(dist) || isstring(dist)
        switch lower(dist)
            case 'unif'
                Z   = GRIDobj(DEM);
                Z.Z = rand(DEM.size);
                if nargin == 3
                    Z.Z = Z.Z .* diff(params);
                    Z.Z = Z.Z + params(1);
                end
            case {'normal','randn'}
                
                Z   = GRIDobj(DEM);
                Z.Z = randn(DEM.size);
                if nargin == 3
                    Z.Z = Z.Z * params(2) + params(1);
                end
            case 'randperm'
                Z   = GRIDobj(DEM);
                Z.Z = reshape(randperm(prod(DEM.size)),DEM.size(1),DEM.size(2));
            case 'randi'
                if nargin == 2
                    params = [1 prod(DEM.size)];
                end
                Z   = GRIDobj(DEM);
                Z.Z  = randi(diff(params),DEM.size) + params(1);
        end
    else
        Z = GRIDobj(DEM);
        Z.Z = random(dist,DEM.size);
    end
end



