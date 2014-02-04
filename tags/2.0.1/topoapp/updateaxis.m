function updateaxis(hax)

hxmin = findobj('Tag','xmin');
hxmax = findobj('Tag','xmax');
hymin = findobj('Tag','ymin');
hymax = findobj('Tag','ymax');
xmin = str2double(get(hxmin,'String'));
xmax = str2double(get(hxmax,'String'));
ymin = str2double(get(hymin,'String'));
ymax = str2double(get(hymax,'String'));
xlims = get(hax,'Xlim');
ylims = get(hax,'Ylim');
if ~isnan(xmin); xlims(1) = xmin;
else set(hxmin,'String','min'); end
if ~isnan(xmax); xlims(2) = xmax;
else set(hxmax,'String','max'); end
if ~isnan(ymin); ylims(1) = ymin;
else set(hymin,'String','min'); end
if ~isnan(ymax); ylims(2) = ymax;
else set(hymax,'String','max'); end
set(hax,'Xlim',xlims,'Ylim',ylims);

end %