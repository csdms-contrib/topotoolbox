function straightcarving(app,IX)
% Straight carving
[p1(1),p1(2)] = ind2sub(app.DEM.size,IX(1));
[p2(1),p2(2)] = ind2sub(app.DEM.size,IX(2));

p = p1;
d = p2-p1;
N = max(abs(d));
s = d/N;

subs = zeros(N,2);
subs(1,:) = p1;

for ii=2:N+1
   p = p+s;
   subs(ii,:) = round(p);
end

IX = sub2ind(app.DEM.size,subs(:,1),subs(:,2));

switch find(cell2mat(get(app.carvemeth,'Value')));
    case 1
        % compute elevation values linearly interpolated
        d  = hypot(subs(:,1)-subs(1,1),subs(:,2)-subs(1,2));
        dt = cumsum(d);
        s  = (app.DEM.Z(IX(end)) - app.DEM.Z(IX(1)))./dt(end);  
        app.DEM.Z(IX) = app.DEM.Z(IX(1)) + s*dt;
        
    case 2
        for r = 2:numel(IX);
            app.DEM.Z(IX(r)) = min(app.DEM.Z(IX(r)),app.DEM.Z(IX(r-1)));
        end
        
end

end