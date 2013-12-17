classdef WATERSHEDobj
    
%    [DB] = WATERSHEDobj(FD,ix)
%    [DB] = WATERSHEDobj(FD,x,y)
    
    
    properties
        outlet
        ix
        x
        y
        distance
    end
    
    
    methods
        
        function [DB] = WATERSHEDobj(FD,varargin)
            narginchk(2,3)
            if nargin==2 % call with index
                DB.outlet = double(varargin{1});
                [y0,x0] = ind2sub(FD.size,DB.outlet);
            else % call with subsrcipts
                x0 = double(varargin{1});
                y0 = double(varargin{2});
                DB.outlet = sub2ind(FD.size,y0,x0);
            end
            
            % extract drainge basin
            db = drainagebasins(FD,DB.outlet);
            DB.ix = find(db.Z);
            XY = bwtraceboundary(db.Z,[y0,x0],'N',8);
            rows = XY(:,1);
            cols = XY(:,2);
            
            % get coordinate pairs
            %[DB.x,DB.y] = sub2coord(FD,DB.ix);
            xy =  [double(rows(:)) double(cols(:)) ones(numel(rows),1)] * FD.refmat;
            DB.x = xy(:,1);
            DB.y = xy(:,2);
                
            DB.distance = getdistance(DB.x,DB.y);
            
        end
        
    end
    
end