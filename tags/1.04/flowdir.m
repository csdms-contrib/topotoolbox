function [M,W0] = flowdir(X,Y,dem,varargin)

% multiple and single flow direction algorithm for Digital Elevation Models
%
% Syntax
%
%     M = flowdir(X,Y,dem)
%     M = flowdir(X,Y,dem,'propertyname',propertyvalue,...)
%     [M,W] = flowdir(...)
%
% Description
% 
%     Multiple and single flowdirection algorithm that routes 
%     through flat terrain (not sinks). Note that this function works
%     best for small grids (smaller 600x600).
%
% Input
%
%     X,Y       coordinate matrices created by meshgrid
%     dem       digital elevation model same size as X and Y
%
% Properties
%
% propertyname     propertyvalues
%
%     'type'            'multi' (default): multiple flowdirection (dinf)
%                       'single': single flow direction (d8). Flow occurs 
%                       only along the steepest descent
% 
%     'exponent'        exponent governing the relation between flow
%                       direction and slope. Default is 1.1, which means,  
%                       there is a nonlinear relation. You may want to 
%                       increase the exponent when flow direction should 
%                       rather follow a steepest descent (single) flow 
%                       direction (e.g. 5). This option is only effective 
%                       for multiple flowdirection.
%
%     'fillsinks'       false (default): pits are not removed
%                       true: pits are filled (requires the Image
%                       Processing Toolbox)
% 
%     'mode'            'default': deterministic flow
%                       'random': totally random flow to downward neighbors
%                       'randomized': deterministic flow with noise
% 
%     'routeflats'      'route' (default), 'cross' or 'none' decides upon 
%                       the method applied to route over flats/plateaus.
%                       'route' uses the function routeflats
%                       'cross' uses the function crossflats
%                       'none' does not apply any routing through flats
% 
%     'edges'           decide on how to handle flow on grid edges. 
%                       'closed' (default) forces all water to remain on 
%                       the grid, 'open' assumes that edge cells loose the  
%                       ratio
%                       r = # of neighbor cells/8
%                       of water.
%
% Output
%
%     M         flowdirection (sparse matrix)
%     W0        amount of "water" in each cell. Usually W is
%               ones(numel(dem),1). Yet, if you use the routeflats option
%               'cross' W includes the area of the flat terrain. You can 
%               include W as an additional parameter to flowacc then.
%
% Example
% 
%     load exampleDEM
%     M = flowdir(X,Y,dem,...
%                   'exponent',5,...
%                   'edges','open',...
%                   'mode','randomized');
% 
%
% Required m-files
% ixneighbors 
%
% See also: EZFLOWACC, FLOWDIR_SINGLE, FLOWACC_LM
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009



if nargin<3;
    error('TopoToolbox:incorrectinput',...
          'wrong number of input arguments')
else
    siz = size(dem);
    if any(size(X) ~= size(Y)) || any((size(X) ~= siz));
        error('TopoToolbox:incorrectinput',...
          'X, Y and dem must have same size')
    end
end

% general values
nrc = numel(dem);

% check input using PARSEARGS
params.type        = {'multi','single'};
params.mode        = {'default','randomized','random'};
params.W0          = ones(siz);
params.fillsinks   = false;
params.exponent    = 1.1;
params.routeflats  = {'route','cross','none'};
params.edges       = {'closed','open'};
params = parseargs(params,varargin{:});


% fill sinks if desired
if params.fillsinks
        dem = fillsinks(dem);
end

% *********************************************************************
% normal multiple FLOW DIRECTION calculation

% calculate maximum slope and slope direction
% find neighbors of cells
[ic1,icd1] = ixneighbors(dem);
e = (dem(ic1)-dem(icd1))./hypot(X(ic1)-X(icd1),Y(ic1)-Y(icd1));


switch params.edges
    case 'open';   
        edgecorrection = histc(ic1,(1:nrc)')./8;
end

e = max(e,0);

% *********************************************************************
% flow direction matrix
M = sparse(ic1,icd1,e,nrc,nrc);


% *********************************************************************
% routing through flats
%
% A flat or plateau exists when two or more adjacent cells have the same 
% elevation. Up to now flow direction indicates for these cells
% no exchange of water.
% The subsequent code first identifies flats and then iteratively removes
% them.

switch params.routeflats
    case 'route'
        [icf,icn] = routeflats(dem,params.type);
        M = sparse(icf,icn,1,nrc,nrc)+M;        
    case 'cross'
        if nargout == 1;
            [icf,icn] = crossflats(dem,params.type);
        else
            [icf,icn,W0] = crossflats(dem,params.type);
        end
        M = sparse(icf,icn,1,nrc,nrc)+M;  
    case 'none'
end


% ******************************************************************
% Randomization of amount transferred to another cell
switch params.mode;
    case 'random'
        M = abs(sprandn(M));
    case 'randomized'
        % randomize coefficient. The higher, the more random
        rc = 0.01;
        M = M + (rc * abs(sprandn(M)));
    otherwise
end
% ******************************************************************
% single flow direction, flow concentration
switch params.type
    case 'single'
        [m,IX2] = max(M,[],2);
        i = m==0;
        IX1 = (1:nrc)';
        IX1(i) = [];
        IX2(i) = [];
        M = sparse(IX1,IX2,1,nrc,nrc);
    otherwise
        if params.exponent ~= 1;
            M = M.^params.exponent;
        end
end

% ******************************************************************
% Row standardization of M only necessary when multiple flow dir
switch params.type
    case 'multi'
        M = spdiags(spfun(@(x) 1./x,sum(M,2)),0,nrc,nrc) * M;
end

% ******************************************************************
% edge correction
switch params.edges
    case 'open';
        M = spdiags(edgecorrection,0,nrc,nrc)*M;
end

% check if second output argument is wanted
if nargout == 2 && ~exist('W0','var')
    W0 = ones(size(dem));
end



% and this is it...



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

