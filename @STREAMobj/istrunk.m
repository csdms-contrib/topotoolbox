function tf = istrunk(S)

%ISTRUNK Determines whether STREAMobj consists of trunk streams
%
% Syntax
%
%     tf = istrunk(S)
%
% Description
%
%     istrunk determines whether a STREAMobj S consists of one or several
%     trunk streams. This is the case if there are as many outlet points as
%     there are channelheads.
%
% Input arguments
%
%     S     STREAMobj
% 
% Output arguments
%
%     tf    true or false
%
% See also: STREAMobj/streampoi
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. February, 2022

tf = numel(streampoi(S,'channelhead','ix')) == ...
     numel(streampoi(S,'outlets','ix'));