function [X,ix,txt,txtname] = GRIDobj2pm(varargin)

%GRIDobj2pm combine several GRIDobj into a predictor matrix
%
% Syntax
%
%     [X,ix,txt,txtname] = GRIDobj2pm(A,B,C,...)
%     [t,ix,txt,txtname] = GRIDobj2pm(A,B,C,...,'ds')
%
% Description
%
%     GRIDobj2pm transforms several instances of GRIDobj to a predictor
%     matrix X. The function removes rows with nans. X(:,m) = A.Z(ix) 
%     where A is the m-th GRIDobj in the input argument list and ix
%     is the linear index that refers to non-nan values in A.Z.
%
%     Note that the predictor matrix does not include an offset, i.e. a
%     column with ones.
%
% Input arguments
%
%     A,B,C,...   spatially aligned GRIDobj with numeric values
%
% Output arguments
%
%     X           predictor matrix
%     ix          index into GRIDobj
%     txt         cell array with variable names of input arguments
%     txtname     cell array with strings obtained from A.name, B.name, ...
%
% Example: A simple landform classification using kmeans
%
%     DEM  = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD   = FLOWobj(DEM,'preprocess','carve');
%     G    = gradient8(DEM);
%     logA = log(flowacc(FD));
%     C    = curvature(DEM);
%     [X,ix,txt,txtname] = GRIDobj2pm(G,C,logA);
%     X = zscore(X);
%     IDX = kmeans(X,5);
%     CL = GRIDobj(DEM);
%     CL.Z(ix) = IDX;
%     imageschs(DEM,CL)
%     
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017


if ischar(varargin{end});
    usetall = true;
    nrgrids = nargin - 1;
else
    usetall = false;
    X   = zeros(numel(varargin{1}.Z),nargin);
    nrgrids = nargin;
end

txt = cell(1,nrgrids);
txtname = cell(1,nrgrids);
for r=1:nrgrids
    
    if r>1
        validatealignment(varargin{r},varargin{1});
    end
    
    if ~usetall
        X(:,r) = double(varargin{r}.Z(:));
    else
        if r == 1
            X = tall(double(varargin{r}.Z(:)));
        else
            X = [X, tall(double(varargin{r}.Z(:)))];
        end
    end
    txt{r} = inputname(r);
    txtname{r} = varargin{r}.name;
end

I  = ~any(isnan(X),2);
ix = (1:numel(varargin{1}.Z))';

X  = X(I,:);
if ~usetall
    ix = ix(I);
else
    ix = ix(gather(I));
end
