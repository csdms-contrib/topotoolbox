function [A,M] = ezflowacc(X,Y,dem,varargin)

% easy to use flow accumulation algorithm for Digital Elevation Models
%
% Syntax
%
%     [A,M] = ezflowacc(X,Y,dem)
%     [A,M] = ezflowacc(X,Y,dem,'propertyname',propertyvalue,...)
%
% Description
% 
%     Easy to use single and multiple flowdirection and flowaccumulation  
%     algorithm that routes through flat terrain (not sinks).
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
%                       direction and slope. Default is 1, which means,  
%                       there is a linear relation. You may want to 
%                       increase the exponent when flow direction should 
%                       rather follow a steepest descent (single) flow 
%                       direction (e.g. 5). This option is only effective 
%                       for multiple flowdirection.
%
%     'fillsinks'       false (default): pits are not removed
%                       true: pits are filled 
% 
%     'mode'            'default': deterministic flow
%                       'random': totally random flow to downward neighbors
%                       'randomized': deterministic flow with noise
% 
%     'W0'              W0 is an initiation grid same size as dem and  
%                       refers to the water in each cell before routing  
%                       through the catchment. By default W0 is 
%                       ones(size(dem)).
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
%     A         flowaccumulation (upslope area) grid
%     M         flowdirection (sparse matrix)
%
%
% Example
% 
%     load exampleDEM
%     A = ezflowacc(X,Y,dem,...
%                   'exponent',5,...
%                   'edges','open',...
%                   'mode','randomized');
%     imagesc(X(1,:),Y(:,2),A); axis image; axis xy
% 
%
% Required m-files
% ixneighbors 
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 12. Januar, 2010


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change threshold of nr of elements
% in the dem. When exceeded flowacc_lm is used
maxnrc = .6*1e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

nrc = numel(dem);

% check input using PARSEARGS
params.type        = {'multi','single'};
params.mode        = {'default','randomized','random'};
params.W0          = ones(siz);
params.fillsinks   = false;
params.exponent    = 1;
params.routeflats  = {'route','cross','none'};
params.edges       = {'closed','open'};
params.waitmessage = {'off','on'};
params = parseargs(params,varargin{:});

if params.fillsinks
        dem = fillsinks(dem);
end

if nrc>maxnrc && ...
   (strncmp(params.type,'multi',1) || ...
    ~strncmp(params.mode,'default',1));    
    A = flowacc_lm(dem,params.W0,...
                       'type',params.type,...
                       'mode',params.mode,...
                       'routeflats',params.routeflats,...
                       'waitmessage',params.waitmessage,...
                       'exponent',params.exponent);

        if nargout == 2;
            warning('TopoToolbox:arraytoolarge',...
              ['could not return multiple flowdirection matrix \n'...
               'since it requires too much memory. Use single \n' ...
               'flow direction algorithm instead.'])
            M = [];
        end
        
elseif nrc>maxnrc && strncmp(params.type,'single',1);
    M  = flowdir_single(dem,'routeflats',params.routeflats,...
                            'edges',params.edges);
    A  = flowacc(M,params.W0);
elseif nrc<=maxnrc;
    M  = flowdir(X,Y,dem,...
                       'type',params.type,...
                       'mode',params.mode,...
                       'routeflats',params.routeflats,...
                       'edges',params.edges,...
                       'exponent',params.exponent);
    A  = flowacc(M,params.W0);
end
    
            
        



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

