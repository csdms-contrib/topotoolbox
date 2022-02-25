function clr = ttclr(feature)

%TTCLR Colors for common map features
%
% Syntax
%
%     clr = ttclr(feature)
%
% Description
%
%     ttclr returns a number of rgb values for common usage. 
%
% Input arguments
%
%     feature     ttclr without input arguments for a list of possible
%                 features.
%
% Output arguments
%
%     clr         rgb-triple 
%
% See also: ttscm, ttcmap
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. February, 2022

validfeatures = {'lake','lakeoutline','river','glacier','desert'};

if nargin == 0
    clr = validfeatures;
    return
end

feature = validatestring(feature,validfeatures,1);

switch lower(feature)
    case 'lake'
        clr = [165, 191, 221]/255;
    case {'lakeoutline','river'}
        clr = [0 0.4470 0.7410];
    case 'glacier'
        % https://support.esri.com/en/technical-article/000010027
        % 10% Gray
        clr = [225 225 225]/255;
    case 'desert'
        % https://support.esri.com/en/technical-article/000010027
        % Sahara sand
        clr = [255 235 190]/255;
end