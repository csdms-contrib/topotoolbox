function setextent(extent,ax)

% set current axis extent 
%
% Syntax
%
%     extent = setextent(extent,ax)
%
% Description
%
%     getextent and setextent are two small wrapper functions (around set
%     and get) to quickly zoom to a specified extent in an axis object.
%
% Example
% 
%     load exampledem
%     imageschs(X,Y,dem)
%     % execute til here and zoom to a desired level using the zoom tool
%     e = getextent;
%     imageschs(X,Y,dem,gradient8(dem))
%     setextent(e)
%
% See also: IMAGESCHS, GETEXTENT
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 11. April, 2011



if nargin == 1
    ax = gca;
end

zoom reset
set(ax,{'xlim','ylim'},extent)