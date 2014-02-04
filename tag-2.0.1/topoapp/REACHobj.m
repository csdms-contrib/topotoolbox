classdef REACHobj
    
    
    properties
        ix
        x
        y
        distance
    end
    
    
    
    methods
        
        function R = REACHobj(FD,varargin)
            
            narginchk(1,3)
            if nargin==2 % call with index
                IX = double(varargin{1});
                ix1 = IX(1);
                ix2 = IX(2);
            else % call with subsrcipts
                X = double(varargin{1});
                Y = double(varargin{2});
                IX = sub2ind(FD.size,Y,X);
                ix1 = IX(1);
                ix2 = IX(2);
            end
            
            [r1,~] = flowpathextract(FD,ix1);
            [r2,~] = flowpathextract(FD,ix2);
            
            % Find overlapping parts
            [c,ia,ib] = intersect(r1,r2,'stable');

            if isempty(c)
                warning('Topoapp:badreach','Found no channelized connection between points.')
            else
                % Assemble reach
                r1(ia) = []; r2(ib) = [];
                R.ix = double([r1;c(1);flipud(r2)]);
                
                % get coordinate pairs
                %[R.x,R.y] = ind2coord(FD,R.ix);
                [rows,cols] = ind2sub(FD.size,R.ix);
                xy =  [double(rows(:)) double(cols(:)) ones(numel(rows),1)] * FD.refmat;
                R.x = xy(:,1);
                R.y = xy(:,2);
                
                R.distance = getdistance(R.x,R.y);
            end
        end
        
%         function h = plot(R,varargin)
% 
%                 h = plot(R.x,R.y,varargin{:});
% 
%         end
        
    end
    
    
    
end