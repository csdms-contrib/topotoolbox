function updateimage(app)
% update display
%
% updateimage(app)
%

set(app.im,'CData',app.DEM.Z);
if any(get(app.im,'AlphaData'))
    set(app.im,'AlphaData',1-app.SINKS.Z*.5);
end

end