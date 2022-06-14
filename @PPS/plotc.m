function [hl,hp] = plotc(P,cova,varargin)

%PLOTC plot a colored stream network
%
% Syntax
%      
%     plotc(P,a)
%     plotc(P,nal)
%     [hl,hp] = plotc(...)
%
% Description
%
%     plotc plots the planform stream network with additional coloring
%     obtained from a GRIDobj or a node attribute list.
%
% Input arguments
%
%     P      Instance of PPS
%     DEM    GRIDobj
%     nal    node attribute list
%
%     Parameter name value pairs
%
%     'linewidth'  Default = 1.5
%
% Output arguments
%
%     hl      handle to surface object
%     hp      handle to points
%
% Example
%
%
%
% See also: PPS/plotc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. November, 2015

p = inputParser;         
p.FunctionName = 'STREAMobj/plotc';
addParameter(p,'linewidth',1.5);
parse(p,varargin{:});

if nargin == 1
    cova = P.S.distance;
else
    cova    = getcovariate(P,cova);
end

[x,y,c] = STREAMobj2XY(P.S,cova);

z = x*0;
ht = surface([x x],[y y],[z z],[c c],...
        'facecolor','none',...
        'edgecolor','flat',...
        'linewidth',p.Results.linewidth);
if nargout >= 1
    hl = ht;
end

if ishold
    keephold = true;
else
    keephold = false;
end

hold on
hp = plotpoints(P);

if ~keephold
    hold off
end

