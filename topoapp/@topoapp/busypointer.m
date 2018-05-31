function busypointer(app,startend)
h = gcf;
if startend
    set(h,'Pointer','watch');
    ch = findobj(h,'Type','uicontrol');
    set(ch,'Enable','off')
    waitfor(ch)
else
    set(h,'Pointer','arrow');
    ch = findobj(h,'Type','uicontrol');
    set(ch,'Enable','on')
    waitfor(ch)
end
end %