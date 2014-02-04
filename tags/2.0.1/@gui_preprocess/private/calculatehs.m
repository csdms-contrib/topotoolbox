function calculatehs(~,~,app)

% fill sinks
app.DEMf = fillsinks(app.DEM);

% recalculate RGB
S = app.SINKDEPTH;
S(S==0) = nan;
app.RGB  = imageschs(app.DEM,S);
if ndims(get(app.im,'Cdata')) == 3;
    set(app.im,'Cdata',app.RGB);
end
end