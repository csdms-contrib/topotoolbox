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
%     'K'      erosional efficiency (node-attribute list or GRIDobj), which
%              may vary spatially. If 'K' is supplied, then chitransform
%              returns the time needed for a signal (knickpoint)
%              propagating upstream from the outlet of S. If K has units
%              m^(1-2m)/y, then time will have units of y. Note that
%              calculating the response time requires the assumption that 
%              n = 1.
%     'plot'   0 : no plot (default)
%              1 : chimap
%     'correctcellsize' {true} or false. If true, than the function will
%              calculate areas in unit^2. This is required if the output of 
%              flowacc is used as input. If the units in A are already 
%              m^2, then set correctcellsize to false.
%     'tribsonly' [] (default) or STREAMobj
%              If a STREAMobj St (must be a subset of S) is supplied, then the
%              function calculates the tributaries in S to St and
%              calculates chi only for these tributaries.
%              
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
% Date: 21. March, 2021


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
addParamValue(p,'tribsonly',[])

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
    a = (1./(K)).*(1./a).^p.Results.mn;
end

if ~isempty(p.Results.tribsonly)
    Scopy = S;
    S = modify(S,'tributaryto2',p.Results.tribsonly);
    a = nal2nal(S,Scopy,a);
end

% cumulative trapezoidal integration
c = cumtrapz(S,a);

if ~isempty(p.Results.tribsonly)
    c = nal2nal(Scopy,S,c,0);
    S = Scopy;
end

% plot if required
if p.Results.plot
    plotc(S,c)
end


