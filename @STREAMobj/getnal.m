function varargout = getnal(S,varargin)

% get node attribute list
%
% Syntax
%
%     z = getnal(S,DEM)
%     [z1,z2,...] = getnal(S,DEM1,DEM2,...)
%     nalstruct = getnal(S,DEM1,DEM2,...,'struct')
%
% Description
%
%     getnal returns node attribute lists (nal) extracted from GRIDobjs 
%     aligned with S.
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
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. June, 2014

if ~ischar(varargin{end})
    varargout = cell(1,max(nargout,numel(varargin)));
    for r = 1:max(nargout,numel(varargin))
        validatealignment(S,varargin{r});
        varargout{r} = double(varargin{r}.Z(S.IXgrid));
    end
else
    for r = 1:numel(varargin)-1;
        validatealignment(S,varargin{r});
        OUT.(inputname(r+1)) = double(varargin{r}.Z(S.IXgrid));
    end
    varargout{1} = OUT;
end
