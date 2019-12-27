function setextent(extent,ax)

% set current axis extent 
%
% Syntax
%
%     setextent(extent,ax)
%     setextent(DEM,ax)
%     setextent(
%
% Description
%
%     getextent and setextent are two small wrapper functions (around set
%     and get) to quickly zoom to a specified extent in an axis object.
%
% Input parameters
%
%     
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
if iscell(extent)
    set(ax,{'xlim','ylim'},extent)
elseif isa(extent,'GRIDobj')
    [x,y] = getoutline(extent);
    set(ax,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
elseif isstruct(extent)
    x = [extent.X];
    y = [extent.Y];
    set(ax,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
elseif isa(extent,'STREAMobj')
    v = info(extent,'boundingbox');
    set(ax,'xlim',v(1:2),'ylim',v(3:4));
elseif isnumeric(extent)
    x = extent(:,1);
    y = extent(:,2);
    set(ax,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
    
end
    