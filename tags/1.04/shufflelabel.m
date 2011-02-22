function L = shufflelabel(L,r)

% shufflelabel randomly relabels a label matrix
%
% Syntax
%
%     L = shufflelabel(L)
%     L = shufflelabel(L,reset)
%
% Description
%
%     shufflelabel randomly changes the order of labels in the  
%     label matrix L. Zeros and nans are ignored. L can be every
%     numerical data type, char or a cell array of strings.
%
%     When called with two input arguments, reset is either true or 
%     false. If true, label values in L are reset to range from
%     one to numel(unique(L(:)). Note that this only works for 
%     numeric arrays.
%
%
% Example
% 
%     L = kron([1 2 3;3 2 3],[1 1; 2 2])
% 
% 
%     L =
% 
%          1     1     2     2     3     3
%          2     2     4     4     6     6
%          3     3     2     2     3     3
%          6     6     4     4     6     6
% 
%     Ls = shufflelabel(L)
% 
%     Ls =
% 
%          3     3     6     6     1     1
%          6     6     2     2     4     4
%          1     1     6     6     1     1
%          4     4     2     2     4     4
% 
%
%     Following expressions work equally well
%     Ls = shufflelabel(uint8(L))
%     Ls = shufflelabel(char(L+80))
%     Ls = shufflelabel(mat2cell(char(L+80),ones(4,1),ones(6,1)))
%
%
% See also: RANDPERM, BWLABEL, LABELMATRIX
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 16. September, 2009



% essentially, error checking is not necessary
error(nargchk(1,2,nargin));

% is L a numeric array? 
inum = isnumeric(L);

% reset labeling?
if nargin == 1;
    % by default no.
    r = false;
else
    % applies only if L is a numeric array.
    r = r&&inum;
end

% size of L
siz = size(L);

% force column vector
L   = L(:);

% exclude zeros and nans from shuffling (only when L is 
% a numeric array)
if inum
    I   = ~(L==0 | isnan(L));
else
    I   = ':';
end

% find unique elements in label vector
[uniqueL,ix,ix] = unique(L(I));

% shuffle labels
if r
    uniqueLS = randperm(numel(uniqueL));
else
    uniqueLS = uniqueL(randperm(numel(uniqueL)));
end

% and map labels back into L
L(I) = uniqueLS(ix);

% finally reshape back to original size
L = reshape(L,siz);

