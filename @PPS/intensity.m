function [in,inci] = intensity(P,alpha)

%INTENSITY calculate intensity (density) of points on the stream network
%
% Syntax
%
%     in = intensity(P)
%     [in,inci] = intensity(P,alpha);
%
% Description
%
%     intensity measures the expected number of points per unit length.
% 
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%     alpha  significance level (default = 0.05)
%
% Output arguments
%
%     in     intensity (scalar)
%     inci   confidence bounds of intensity (two element vector)
%
% See also: PPS, PPS/npoints 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

n  = npoints(P);
l  = tlength(P);
in = n/l;   
if nargout > 1
    % If error bounds required
    if nargin == 1
        alpha = 0.05;
    else
        validateattributes(alpha,{'numeric'},{'>',0,'<',1},...
            'intensity','alpha',2);
    end
    % Fit loglinear model with intercept only
    mm = fitloglinear(P,getnal(P.S)+1,'intercept',false);
    [p,pci]   = predict(mm,1,'alpha',alpha);
    d   = l./numel(P.S.x);
    in = p./d; %.cellsize;   
    inci = pci./d;
end