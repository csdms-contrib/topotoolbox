classdef PROFILEobj
    
%    [P] = PROFILEobj(DEM,ix)
%    [P] = PROFILEobj(DEM,x,y)
    
    
    properties
        name
        x
        y
        distance
    end
    
    
    methods
        
        function [P] = PROFILEobj(DEM,x,y)
            [P.x,P.y] = sub2coord(DEM,y,x);
            [P.distance] = getdistance(P.x,P.y);
        end
        
%         function h = plot(P,varargin)
%             % plot object trace in map
% %             hold on
% %             if ~isempty(varargin)
%                 h = plot(P.x,P.y,varargin{:});
% %             else
% %                 h = plot(P.x,P.y,'color',P.color);
% %             end
% %             hold off
% %             set(h,'Tag',class(P))
% %             if ~isempty(P.name)
% %                 set(h,'DisplayName',P.name)
% %             end
% %             P.handle = h;
% %             P.visible = 1;
%         end
        
    end
    
end