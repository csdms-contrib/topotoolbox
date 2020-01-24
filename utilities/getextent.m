function extent = getextent(ax)

%GETEXTENT get current axis extent 
%
% Syntax
%
%     extent = getextent(ax)
%
% Description
%
%     getextent and setextent are two small wrapper functions (around set
%     and get) to quickly zoom to a specified extent in an axis object.
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM);
%     % execute til here and zoom to a desired level using the zoom tool
%     e = getextent;
%     imageschs(DEM,gradient8(DEM))
%     setextent(e)
%
% See also: IMAGESCHS, SETEXTENT, PADEXTENT
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 23. January, 2020

if nargin == 0
    ax = gca;
end
extent = get(ax,{'xlim','ylim'});  % Get axes limits.

