function extent = getextent(ax)

% get current axis extent 
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
%     load exampledem
%     imageschs(X,Y,dem);
%     % execute til here and zoom to a desired level using the zoom tool
%     e = getextent;
%     imageschs(X,Y,dem,gradient8(dem))
%     setextent(e)
%
% See also: IMAGESCHS, SETEXTENT
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 11. April, 2011

if nargin == 0
    ax = gca;
end
extent = get(ax,{'xlim','ylim'});  % Get axes limits.

