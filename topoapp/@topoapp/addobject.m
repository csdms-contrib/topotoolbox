function [app] = addobject(app,objdata,varargin)
% ADDOBJECT adds object data to instance of class TOPOAPP
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

% Parse inputs
p = inputParser;
p.FunctionName = 'addobject';
addRequired(p,'app',@(x) isa(x,'topoapp'));
addParameter(p,'name','',@(x) ischar(x))
addParameter(p,'visible',false,@(x) islogical(x))
addParameter(p,'handle',[]);%,@(x) isnumeric(x))
parse(p,app,varargin{:});

objectname = class(objdata);
handle = p.Results.handle;
name = p.Results.name;
visible = p.Results.visible;
ct = app.objects.(objectname).ct + 1;
if isempty(name)
    name = [objectname,' ',num2str(ct)];
end
ix = length(app.objects.(objectname).data)+1;
app.objects.(objectname).names{ix} = name;
app.objects.(objectname).handles(ix) = handle;
app.objects.(objectname).visible(ix) = visible;
app.objects.(objectname).data{ix} = objdata;
app.objects.(objectname).ct = ct;
end %