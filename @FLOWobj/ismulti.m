function tf = ismulti(FD,typetest)
% check if FD is multi or single flow direction
%
% Syntax
%
%     tf = ismulti(FD)
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016

if nargin == 1
    typetest = true;
end

if typetest
tf = strcmp(FD.type,'multi');
else 
tf = any(histcounts(FD.ix,1:(prod(FD.size)+1))>1);
end
