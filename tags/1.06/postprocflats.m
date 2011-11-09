function A = postprocflats(flats,A,fun)

% postprocess flat terrain for visualization purpose
%
% Syntax
%
%     Ap = postprocflats(flats,A,fun)
%
% Description
%
%     Sometimes it may be nice for visualizations to manipulate the flow 
%     accumulation raster so that flat areas such as lakes have constant
%     values instead drainage patterns returned by flow direction.
%     postprocsinks assigns equal values to each connected, flat area using
%     a user-defined function.
%
% Input
%
%     flats     logical array indicating flats (as returned by the function
%               identifyflats)
%     A         flow accumulation array
%     fun       scalar, string or function handle of a function that takes 
%               a vector and returns a scalar such as min, max, median, .
%               mean etc (default = @max). To assign the same value to each
%               flat, use following syntax: A2 = postprocflats(I,A,100000);
% 
% Output
%
%     Ap        postprocessed flow accumulation array
%
% Example
%
%     [X,Y,dem] = peaks(200);
%     dem       = fillsinks(dem);
%     A         = ezflowacc(X,Y,dem,'type','single');
%     I         = identifyflats(dem);
%     A         = postprocflats(I,A,'max');
%     surf(X,Y,dem,log(A));
%
%
% See also: IDENTIFYFLATS, FUNCTION_HANDLE
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 1. June, 2010


% check input arguments
if ~isa(flats, 'logical')
    flats = flats>0 & ~isnan(flats);
end

if ~isequal(size(flats),size(A))
    error('TopoToolbox:input','flats and A must have same size')
end

if nargin==3
    
    if isa(fun, 'char');
        fun   = str2func(['@(x)' fun '(x)']);
    end
    if isa(fun, 'numeric');
        fun   = str2func(['@(x)' num2str(fun)]);
    end
    
else
    fun   = @max;
end

% find connected components
CC    = bwconncomp(flats,8);
% extract values in A
cA  = cellfun(@(x) A(x),CC.PixelIdxList,'UniformOutput',false);

% compute new values
try
    val = cell2mat(cellfun(fun,cA,'UniformOutput',false));
catch %#ok
    error('the function must take a vector and return a scalar')
end


% and write them in A
for r = 1:numel(CC.PixelIdxList);
    A(CC.PixelIdxList{r}) = val(r);
end