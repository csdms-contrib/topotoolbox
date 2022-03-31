function tf =  isempty(S)

%ISEMPTY Determine whether a STREAMobj is empty
%
% Syntax
%
%     tf = isempty(S)
%
% Description
%
%     isempty determines if S has at least two nodes and one edge. 
%
% Input arguments
%
%     S   STREAMobj
%
% Output arguments
%
%     tf  true or false
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 31. March, 2022

tf = isempty(S.x);