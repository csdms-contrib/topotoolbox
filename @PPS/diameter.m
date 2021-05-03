function [dmax,IX] = diameter(P,varargin)

%DIAMETER Returns the maximum possible distance in a stream network
%
% Syntax
%
%      dmax = diameter(P)
%      dmax = diameter(p,'val',d)
%
% Description
%
%      DIAMETER returns the diameter of the stream network. In graph
%      theory, the diameter of a network refers to the maximum distance
%      between any two pairs of nodes. 
%
% Input arguments
%
%      P     Instance of PPS
%
%      Parameter name/value pairs
%
%      'd3d'   {false} or true: distances in 3d
%      'val'   node attribute list with custom distance value. By default,
%              the distance is measured in mapunits
%      'usepoints' {false} or true: If false, the diameter is calculated
%              for the entire network, if true, only the distance between 
%              the points is calculated.
%      'plot'  {false} or true: If true, the function will create a plot
%              with the start and end point of the path with the maximum 
%              length.
%
% Output arguments
%
%      dmax    maximum distance 
%      IX      linear index into the GRIDobj with the end points of the
%              path with the maximum length.
% Example
%
%
%
%
% See also: PPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 22. December, 2020


p  = inputParser;
addParameter(p,'d3d',false);
addParameter(p,'val',[]);
addParameter(p,'plot',false);
addParameter(p,'usepoints',false);
parse(p,varargin{:});

Results = p.Results;
plotit  = Results.plot;
usepoints = Results.usepoints;
Results = rmfield(Results,{'plot' 'usepoints'});

if ~usepoints
    P = generatepoints(P,'type',{'outlet','channelhead'});
end
d  = pointdistances(P,Results);
d(isinf(d)) = -inf;
[dmax,ix] = max(d(:));
[r,c]  = ind2sub(size(d),ix);
P.PP = P.PP([r c]);
P.PP = P.PP(:);
IX    = P.S.IXgrid(P.PP);

if plotit
    plot(P.S,'k');
    % get shortest path
    path = shortestpath(as(P,'graph'),P.PP(1),P.PP(2));
    hold on
    plot(P.S.x(path),P.S.y(path),'r');
    plotpoints(P);
    hold off
end
