function tf = ismulti(M)

% determine whether matrix is a multiple flow direction matrix 
%
% Syntax
% 
%     tf = ismulti(M)
%
% Description
%
%     Various functions in the TopoToolbox only work with single flow
%     direction matrices. ismulti checks if the flow direction matrix is a
%     multiple flow direction matrix.
%
% Input
%
%     M         flow direction matrix
%
% Output
%
%     tf        true if M is a multiple flow direction matrix
%
% See also: multi2single, flowdir
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 16. April, 2010


tf = any(full(sum(spones(M)~=0,2)>1));