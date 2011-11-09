function [M,W0] = flowdir_single(DEM,varargin)

% single flow direction (upslope area)
% 
% Syntax
%
%     M = flowdir_single(dem)
%     M = flowdir_single(dem,'propertyname',propertyvalue,...)
%     [M,W] = flowdir_single(...)
%
% Description
%
%     flowdir_single calculates the flow direction matrix M of a digital 
%     elevation model. The algorithm employs the single flowdirection 
%     algorithm that routes water down the steepest descend. Compared to 
%     ezflowacc and flowdir, flowdir_single is memory efficient and can 
%     handle large digital elevation models.
%
% Input
%
%     dem       Digital Elevation Model
%
% propertyname     propertyvalues
%
%     'routeflats'      choose method to route over flats/plateaus.
%                       'route' uses the function routeflats
%                       'cross' uses the function crossflats
%                       'geodesic' uses routegeodesic (requires Matlab
%                       2011b or higher)
%                       'none' does not apply any routing through flats
%     'edges'           decide on how to handle flow on grid edges. 
%                       'closed' (default) forces all water to remain on 
%                       the grid, 'open' assumes that all water on edges
%                       leaves the grid.
%
% 
% Output
%
%     M         sparse flowdirection matrix
%     W0        amount of "water" in each cell. Usually W is
%               ones(numel(dem),1). Yet, if you use the routeflats option
%               'cross' W includes the area of the flat terrain. You can 
%               include W as an additional parameter to flowacc then.
%               
%
% Example
% 
%     load exampleDEM
%     M = flowdir_single(dem);
%     gplot(M,[X(:) Y(:)]);
%     axis image
%     title('flow network')
%
% See also: FLOWDIR, EZFLOWACC, FLOWACC_LM
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009




% General variables
% cellsize
d      = 1; 

siz    = size(DEM);
nrc    = numel(DEM);
IXc    = reshape((1:nrc)',siz);

% check input using PARSEARGS
params.routeflats  = {'route','cross','geodesic','none'};
params.edges       = {'closed','open'};
params.diagonal    = {'hypot','cs'};
params = parseargs(params,varargin{:});

% diagonal distance
switch lower(params.diagonal)
    case 'hypot'
    dd = hypot(d,d);
    case 'cs'
        dd = d;
end

% lower the outer rim of the dem to ensure flow towards the edges 
B = true(siz);
B(2:end-1,2:end-1) = false;
DEM(B) = DEM(B)-1e-4; % arbitrary!!!
DEM    = padarray(DEM,[1 1],NaN);

% anonymous functions for neighborhood operations
neighfun = cell(8,2);

neighfun{1,1} = @() (DEM(2:end-1,2:end-1)-DEM(1:end-2,2:end-1));
neighfun{2,1} = @() (DEM(2:end-1,2:end-1)-DEM(1:end-2,3:end))/dd;
neighfun{3,1} = @() (DEM(2:end-1,2:end-1)-DEM(2:end-1,3:end));
neighfun{4,1} = @() (DEM(2:end-1,2:end-1)-DEM(3:end,3:end))/dd;

neighfun{5,1} = @() (DEM(2:end-1,2:end-1)-DEM(3:end,2:end-1));
neighfun{6,1} = @() (DEM(2:end-1,2:end-1)-DEM(3:end,1:end-2))/dd;
neighfun{7,1} = @() (DEM(2:end-1,2:end-1)-DEM(2:end-1,1:end-2));
neighfun{8,1} = @() (DEM(2:end-1,2:end-1)-DEM(1:end-2,1:end-2))/dd;

neighfun{1,2} = @() -1;
neighfun{2,2} = @() -1 + siz(1);
neighfun{3,2} = @() siz(1);
neighfun{4,2} = @() 1 + siz(1);

neighfun{5,2} = @() 1;
neighfun{6,2} = @() 1 - siz(1);
neighfun{7,2} = @() - siz(1);
neighfun{8,2} = @() -1 - siz(1);


% search maximum downward neighbor
IXn = zeros(siz);
G   = IXn;

% use randperm to reduce bias due to preferential search direction
% ix = randperm(8);
ix = 1:8;

for neigh = ix;
    G2       = neighfun{neigh,1}();
    I        = G2>G;
    G(I)     = G2(I);
    IXn(I)   = IXc(I)+neighfun{neigh,2}();
end

DEM = DEM(2:end-1,2:end-1);

G   = G(:);
IXn = IXn(:);
IXc = IXc(:);

I = G==0 | isnan(G);
IXc(I) = [];
IXn(I) = [];


% routing through flats
%
% A flat or plateau exists when two or more adjacent cells have the same 
% elevation. Up to now flow direction indicates for these cells
% no exchange of water.
% The subsequent code first identifies flats and then iteratively removes
% them.

switch params.routeflats
    case 'route'
        [icf,icn] = routeflats(DEM,'single');    
    case 'cross'
        if nargout == 1;
            [icf,icn] = crossflats(DEM,'single');
        else
            [icf,icn,W0] = crossflats(DEM,'single');
            W0t = ones(siz);
            W0t(icf) = W0;
            W0  = W0t;
        end
    case 'geodesic'
        [icf,icn] = routegeodesic(DEM,'single');
    case 'none'
end


% handle edges
switch params.edges
    case 'open'
        B = true(siz);
        B(2:end-1,2:end-1) = false;
        I = B(IXc) & B(IXn);
        IXc(I) = [];
        IXn(I) = [];
end

% create sparse flow direction matrix
switch params.routeflats
    case 'none'
        M = sparse(IXc,IXn,1,nrc,nrc);
    otherwise
        switch params.edges
            case 'closed'
                M = sparse([IXc;icf],[IXn;icn],1,nrc,nrc);
            case 'open'
                B = true(siz);
                B(2:end-1,2:end-1) = false;
                IXc = [IXc;icf];
                IXn = [IXn;icn];
                I = B(IXc) & B(IXn);
                IXc(I) = [];
                IXn(I) = [];
                M = sparse(IXc,IXn,1,nrc,nrc);
            otherwise
        end
end

if nargout == 2 && ~exist('W0','var')
     W0 = ones(siz);
end

% ... and that's it



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

function X = parseargs(X,varargin)
%PARSEARGS - Parses name-value pairs
%
% Behaves like setfield, but accepts multiple name-value pairs and provides
% some additional features:
% 1) If any field of X is an cell-array of strings, it can only be set to
%    one of those strings.  If no value is specified for that field, the
%    first string is selected.
% 2) Where the field is not empty, its data type cannot be changed
% 3) Where the field contains a scalar, its size cannot be changed.
%
% X = parseargs(X,name1,value1,name2,value2,...) 
%
% Intended for use as an argument parser for functions which multiple options.
% Example usage:
%
% function my_function(varargin)
%   X.StartValue = 0;
%   X.StopOnError = false;
%   X.SolverType = {'fixedstep','variablestep'};
%   X.OutputFile = 'out.txt';
%   X = parseargs(X,varargin{:});
%
% Then call (e.g.):
%
% my_function('OutputFile','out2.txt','SolverType','variablestep');

% The various #ok comments below are to stop MLint complaining about
% inefficient usage.  In all cases, the inefficient usage (of error, getfield, 
% setfield and find) is used to ensure compatibility with earlier versions
% of MATLAB.

remaining = nargin-1; % number of arguments other than X
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
% Take input arguments two at a time until we run out.
while remaining>=2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname,fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X,fieldname); %#ok
        newvalue = varargin{count+1};
        if iscell(oldvalue)
            % Cell arrays must contain strings, and the new value must be
            % a string which appears in the list.
            if ~iscellstr(oldvalue)
                error(sprintf('All allowed values for "%s" must be strings',fieldname));  %#ok
            end
            if ~ischar(newvalue)
                error(sprintf('New value for "%s" must be a string',fieldname));  %#ok
            end
            if isempty(find(strcmp(oldvalue,newvalue))) %#ok
                error(sprintf('"%s" is not allowed for field "%s"',newvalue,fieldname));  %#ok
            end
        elseif ~isempty(oldvalue)
            % The caller isn't allowed to change the data type of a non-empty property,
            % and scalars must remain as scalars.
            if ~strcmp(class(oldvalue),class(newvalue))
                error(sprintf('Cannot change class of field "%s" from "%s" to "%s"',...
                    fieldname,class(oldvalue),class(newvalue))); %#ok
            elseif numel(oldvalue)==1 & numel(newvalue)~=1 %#ok
                error(sprintf('New value for "%s" must be a scalar',fieldname));  %#ok
            end
        end
        X = setfield(X,fieldname,newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
% Check that we had a value for every name.
if remaining~=0
    error('Odd number of arguments supplied.  Name-value pairs required');
end

% Now find cell arrays which were not modified by the above process, and select
% the first string.
notmodified = find(~modified);
for i=1:length(notmodified)
    fieldname = fields{notmodified(i)};
    oldvalue = getfield(X,fieldname); %#ok
    if iscell(oldvalue)
        if ~iscellstr(oldvalue)
            error(sprintf('All allowed values for "%s" must be strings',fieldname)); %#ok
        elseif isempty(oldvalue)
            error(sprintf('Empty cell array not allowed for field "%s"',fieldname)); %#ok
        end
        X = setfield(X,fieldname,oldvalue{1}); %#ok
    end
end




