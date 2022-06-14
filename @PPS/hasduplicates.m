function [tf,I] = hasduplicates(P,plotit)

%HASDUPLICATES checks whether there are duplicate points
%
% Syntax
%
%     tf = hasduplicates(P)
%     tf = hasduplicates(P,plotit)
%     [tf,I] = ...
%
% Description
%
%     hasduplicates checks whether the point pattern P has duplicate point
%     entries. That means, two or more points share the same location.
%
% Input arguments
%
%     P        point pattern on stream network (class PPS)
%     plotit   {false} or true. If true, hasduplicates will plot the
%              network and points and highlight the duplicate points.
%
% Output arguments
%
%     tf     true or false
%     I      returns a logical vector of size [npoints(P) 1] with elements
%            set to true if they are not unique.
%
%
% See also: PPS, PPS/removeduplicates 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

[u,~,b] = unique(P.PP);

if numel(u) == numel(P.PP)
    tf = false;
    I  = false(size(P.PP));
    return
end

tf = true;
np = npoints(P);
n  = accumarray(b,ones(np,1),[numel(u) 1],@sum);
I  = n(b) > 1;

if nargin > 1
    if plotit
        if ishold
            keephold = true;
        else
            keephold = false;
        end
        plot(P.S);
        hold on
        plotpoints(P,'marks',+I)
        
        if ~keephold
            hold off
        end
    end
end
        