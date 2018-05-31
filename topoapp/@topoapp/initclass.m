function [app] = initclass(app,newclass,color)
% INITCLASS creates a new field in TOPOAPP for objects of NEWCLASS that is
% used to store topoapp object data.
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013
classes = app.objects.classes;
if sum(strcmp(classes,newclass))==0
    app.objects.classes{length(classes)+1} = newclass;
    app.objects.(newclass).ct = 0;
    app.objects.(newclass).color = color;
    app.objects.(newclass).names = [];
    app.objects.(newclass).handles = [];
    app.objects.(newclass).visible = [];
    app.objects.(newclass).data = [];
end
end %