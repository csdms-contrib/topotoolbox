function S = onl2struct(ONL,varargin)
%ONL2STRUCT    convert orderednanlist to structure
%
% Syntax
%
%     S = onl2struct(ONL)
%     S = onl2struct(ONL,'attributename1',attribute1,...
%                        'attributename2',attribute2,...)
%
% Description
%
%     ONL2STRUCT converts an orderednanlist to a structure, in which each
%     entry corresponds to the values found between bounding nans.
%
% Input arguments
%
%     ONL    orderednanlist to structure (nan 1 2 3 nan 4 5 6 nan 3 3 ..)
%
%     varargin: pairs of 'attributename' and attribute, where each
%           attribute needs to have the same length as the ONL and nan
%           at the same places
%
% Output arguments
%
%     S      instance of class structure with fields:
%            .IX - the values in ONL
%            .st - the values at the end of each entry
%            .nr - the number of values in each entry
%            optional:
%            .attributename1 - the values of attribute1
%            .(...)
%
% Examples
%
%
% See also:
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


S = struct;
ix = [0;find(isnan(ONL))];
for i = 1 : length(ix)-1
    S(i).IX = ONL(ix(i)+1:ix(i+1));
    S(i).st = S(i).IX([1,end-1])';
    S(i).nr = numel(S(i).IX)-1;
end

for k = 1 : 2 : length(varargin)
    for i = 1 : length(ix)-1
        S(i).(varargin{k}) = varargin{k+1}(ix(i)+1:ix(i+1));
    end
end

end