function varargout = getnal(S,varargin)

%GETNAL get node-attribute list
%
% Syntax
%
%     z = getnal(S,DEM)
%     [z1,z2,...] = getnal(S,DEM1,DEM2,...)
%     nalstruct = getnal(S,DEM1,DEM2,...,'struct')
%     z = getnal(S)
%
% Description
%
%     getnal returns node-attribute lists (nal) extracted from GRIDobjs 
%     aligned with S. A node-attribute list is a column vector that has as 
%     many elements as nodes in the stream network S. These vectors can
%     store all kinds of information related to each node of the stream
%     network such as elevation, upstream area, and gradient. Numerous
%     STREAMobj methods (e.g. gradient or smooth) accept nals as input.
%
%     nals are intricately linked to the network topology. To export nals
%     to nan-punctuated vectors that can be plotted by custom functions,
%     use the function STREAMobj2XY.
%
%     With only S as input argument, getnal returns a node-attribute list
%     of zeros (double).    
%
% Input
%
%     S          STREAMobj
%     DEM        GRIDobj or several GRIDobjs
%     or 
%     (DEM1,DEM2,...)   
%     
% Output
%
%     z (z1,z2)  vector same size as S.IXgrid that contains the values
%                extracted for each node of the stream network
%     nalstruct  structure array with node attribute lists where fieldnames
%                refer to the workspace variable names.
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A = flowacc(FD);
%     S  = STREAMobj(FD,A>500);
%     S = klargestconncomps(S);
%     nal = getnal(S,DEM,A,'struct');
%     chiplot(S,nal.DEM,nal.A);
%
% See also: STREAMobj/isnal
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017

if nargin == 1
    % return nal of zeros (double)
    varargout{1} = zeros(size(S.x));
    return
end

if ~ischar(varargin{end})
    varargout = cell(1,max(nargout,numel(varargin)));
    for r = 1:max(nargout,numel(varargin))
        validatealignment(S,varargin{r});
        varargout{r} = double(varargin{r}.Z(S.IXgrid));
    end
else
    for r = 1:(numel(varargin)-1)
        validatealignment(S,varargin{r});
        OUT.(inputname(r+1)) = double(varargin{r}.Z(S.IXgrid));
    end
    varargout{1} = OUT;
end
        
    