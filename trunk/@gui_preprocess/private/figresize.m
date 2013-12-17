function figresize(~,~,app)
    % Figure resize callback
    fhpos = get(app.fh,'Position');
    set(app.ax,'Position',[0 0 fhpos(3)-200 fhpos(4)]);
    set(app.mpanel,'Position',[fhpos(3)-200 1 199 fhpos(4)-1]);
end