function c = chitransform(S,A,varargin)

%CHITRANSFORM Coordinate transformation using the integral approach
%
% Syntax
%
%     c = chitransform(S,A)
%     c = chitransform(S,A,pn,pv,...)
%
% Description
%
%     CHITRANSFORM transforms the horizontal spatial coordinates of a river
%     longitudinal profile using an integration in upstream direction of
%     drainage area (chi, see Perron and Royden, 2013).
%
% Input arguments
%
%     S     STREAMobj
%     A     upslope area as returned by the function flowacc
%     
% Parameter name/value pairs
%
%     'a0'     reference area (default=1e6)
%     'mn'     mn-ratio (default=0.45)
%     'plot'   0 : no plot (default)
%              1 : chimap
%
% Output argument
%
%     c     node attribute list (nal) of chi values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     A = flowacc(FD);
%     c = chitransform(S,A,'mn',0.45);
%     plotc(S,c)    
%
% See also: chiplot
%
% References:
%     
%     Perron, J. & Royden, L. (2013): An integral approach to bedrock river 
%     profile analysis. Earth Surface Processes and Landforms, 38, 570-576.
%     [DOI: 10.1002/esp.3302]
%
%     TopoToolbox blog posts >>Chimaps in a few lines of codes<<
%     <a href="https://topotoolbox.wordpress.com/2017/08/18/chimaps-in-a-few-lines-of-code-final/">See overview here.</a>
%     
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 29. December, 2015


% Parse Inputs
p = inputParser;         
p.FunctionName = 'chitransform';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj') || isnal(S,x));

addParamValue(p,'mn',0.45,@(x) isscalar(x) || isempty(x));
addParamValue(p,'a0',1e6,@(x) isscalar(x) && isnumeric(x));
addParamValue(p,'plot',false);
addParamValue(p,'correctcellsize',true,@(x) isscalar(x));
addParamValue(p,'K',[],@(x) isempty(x) || isnal(S,x) || isa(x,'GRIDobj'));

parse(p,S,A,varargin{:});

% get node attribute list with elevation values
if isa(A,'GRIDobj')
    validatealignment(S,A);
    a = getnal(S,A);
elseif isnal(S,A)
    a = A;
else
    error('Imcompatible format of second input argument')
end

if ~isempty(p.Results.K)
    calcwithk = true;
    if isnal(S,p.Results.K)
        K = p.Results.K;
    else 
        K = getnal(S,p.Results.K);
    end
else
    calcwithk = false;
end
        
       
if p.Results.correctcellsize
    a = a.*S.cellsize^2;
end

if ~calcwithk
    a = ((p.Results.a0) ./a).^p.Results.mn;
else
    % This transformation is only possible if we assume that n in the
    % mn-ratio is one.
    a = (1./(K.^p.Results.mn)).*((p.Results.a0)./a).^p.Results.mn;
end
c = cumtrapz(S,a);

if p.Results.plot
    plotc(S,c)
end

