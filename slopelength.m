function [SL,SLF] = slopelength(M,X,Y,dem,varargin)

% slope length of a digital elevation model
%
% Syntax
%
%     SL = slopelength(M,X,Y,dem)
%
% Description
%
%     influencemap returns a logical matrix masking the downslope part
%     of the digital elevation model that is drained by specified cells.
%
% Input
%
%     M         single flow direction matrix
%     X,Y       coordinate matrices
%     dem       digital elevation model
%
% Output
%
%     SL        slope length raster  
%
% Example
%
%     load exampleDEM
%     M = flowdir_single(dem);
%     SL = slopelength(M,X,Y,dem)
%     imageschs(X,Y,dem,SL)
%
%
% See also: FLOWDIR
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 4. August, 2010




% error checking
M = spones(M);
if any(sum(M,2)>1);
    error('TopoToolbox:incorrectinput',...
          'single flow direction must be supplied')
end


% defaults
params.distance = {'2D','3D'};
params.type = 'McCool';
params = parseargs(params,varargin{:});


cs  = abs(Y(1)-Y(2));
siz = size(dem);

% calculate upstream flowdistance
switch lower(params.distance)
    case '2d'
        D = flowdistance(M,X,Y);
    case '3d'
        D = flowdistance(M,X,Y,dem);
end




% Calculation of Slope Length
[ic,icd]   = find(M);
[dummy,ix] = sort(D(ic),'descend');

clear dummy

ic          = ic(ix);
icd         = icd(ix);
[icdd,icdd] = ismember(icd,ic); 


% initialize length-slope
SL  = zeros(siz);
% Find non-receivers (cells on ridges)
NG  = sum(M,1)' == 0;
% slope length in non-receiver cells in half the
% cellsize
SL(NG) = cs/2;

NG  = NG(ic);
NG  = find(NG);
% nr of nonreceivers
nNG = numel(NG);

NGix = 1;

D = abs(D(ic)-D(icd));

IX  = 1;


% not vectorized
while NGix <= nNG;

    % distance
    ln = D(IX) + SL(ic(IX));

    if ln > SL(icd(IX));
        SL(icd(IX)) = ln;
        
        % next row to work on
        IX = icdd(IX);
        
        flaggoon = IX~= 0;
    else
        flaggoon = false;       
    end
    
    if ~flaggoon;
        
        NGix = NGix+1;
        
        if NGix>nNG;
        else
            IX  = NG(NGix);
        end
    end
end




% Compute different slope length factors

if nargout == 2

% we need the gradient
G = gradient8(dem);

switch lower(params.type)
        
    case 'mccool'
        % length-slope factor (McCool)
        I  = G < 0.09;
        G = atan(G);
        % ratio of rill to interrill erosion
        f = sin(G)./(0.0896 * (3*sin(G).^0.8 + .56));
        a = 1;
        m = a*f./ (1+a*f);
        
        SLF = zeros(siz);
        
        SLF(I)  = (SL(I)/22.13).^m(I) .* (10.8 * sin(G(I)) + 0.03);
        SLF(~I) = (SL(~I)/22.13).^m(~I) .* (16.8 * sin(G(~I)) - 0.5);
        
    case 'usle'
        % LS factor according to Wischmeier and Smith 1978
        
        G = atan(G);       
        SLF = (SL/22.13).^t * (65.4*sin(G).^2 + 4.56 * sin(G) + 0.0654);
end
 
end

end





%
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
end