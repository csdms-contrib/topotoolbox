function showtools(hObject,~,app)

if app.menu.menucarve == hObject
    set(app.panelfill,'Visible','off')
    set(app.panelcarve,'Visible','on')
    set(app.menu.menucarve,'Checked','on')
    set(app.menu.menufill,'Checked','off')
else
    set(app.panelfill,'Visible','on')
    set(app.panelcarve,'Visible','off')
    set(app.menu.menucarve,'Checked','off')
    set(app.menu.menufill,'Checked','on')
end
end