function [c,l,cc] = histogram(P,varargin)

%HISTOGRAM histogram of point pattern on stream network
%
% Syntax
%
%     c = histogram(P)
%     c = histogram(P,pn,pv,...)
%     [c,label] = ...
%
% Description
%
%     histogram computes the frequency of points on a network. The function
%     subdivides the network into a number of network cells (or bins). The
%     parameter seglength determines the approximate length of each network
%     cell using the function STREAMobj/labelreach. However, only rarely
%     can the network be subdivided into equal-length cells because the
%     functions splits the network at confluences so that individual cells
%     may be shorter than the distance given by seglength. histogram thus
%     does not return the number of points (counts) by default, but the
%     estimate of the probability density function (pdf).
%
% Input arguments
%
%     P      instance of PPS
%
%     Parameter name/value pairs
%
%     'seglength'      cell (or bin) distance. Default is 
%                      tlength(P)/npoints(P)*10
%     'normalization'  'pdf' (default) or 'counts'
%     'voronoi'        false (default) or true. If true, cells are computed
%                      using a voronoi tesselation. 
%
% Output arguments
%
%     c      node-attribute list with pdf/count values
%     label  node-attribute list with labels indicating the cells (bins)
%
% Example
%
%
%
% See also: PPS, PPS/density, PPS/intensity
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

% Check input arguments
p = inputParser;
p.FunctionName = 'PPS/histogram';
addParameter(p,'seglength',tlength(P)/npoints(P)*10);
addParameter(p,'normalization','pdf');
addParameter(p,'voronoi',true);
addParameter(p,'shufflelabel',false);
% Parse
parse(p,varargin{:});

if ~p.Results.voronoi
    l  = labelreach(P.S,'seglength',p.Results.seglength);
else
    d = P.S.distance;
    I = gradient(P.S,mod(d,p.Results.seglength)) < 0 ...
        | streampoi(P.S,'outl','logical');
    PV = PPS(P.S,'PP',I);
    l = voronoi(PV);
end


numl = max(l);
nal = points(P,'nal');

c   = accumarray(l,nal,[numl 1],@sum);

switch lower(p.Results.normalization)
    case 'pdf'  
        if ~p.Results.voronoi
            d   = P.S.distance;
            ld  = accumarray(l,d,[numl 1],@range);
        else
            dd  = zeros(size(d));
            dd(P.S.ix) = d(P.S.ix)-d(P.S.ixc);
            ld  = accumarray(l,dd,[numl 1],@sum);
        end
        c   = c./ld;
end
if nargout == 3
    cc = c;
end

c   = c(l);

if p.Results.shufflelabel
    [uniqueL,~,ix] = unique(l);
    uniqueLS = randperm(numel(uniqueL));
    l = uniqueLS(ix);
end


