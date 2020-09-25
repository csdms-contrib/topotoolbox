function crsapp(S,DEM)

%CRSAPP interactive smoothing of river long profiles
%
% Syntax
%
%     crsapp(S,DEM)
%
% Description
%
%     CRSAPP is an interactive tool to visually assess the results of the
%     function STREAMobj/crs. You can export the results and the parameters
%     to the workspace.
%
%     The graphical user interface allows you to adjust the crs-parameter
%     K, tau, and mingradient with a number of sliders. The K-slider
%     adjusts the degree of smoothing and uses a logarithmic scaling.
%     Common values range between 1 and 10. mingradient adjust the minimum
%     gradient that a profile must go downward in downstream directions.
%     Choose values so that they do not exceed the true downstream
%     gradients. tau-values range between 0.0001 and 0.9999 and refer to
%     the quantile along which the smoothed profile should run along the
%     measured profile.
%
% Input parameters
%
%     S      STREAMobj
%     DEM    digital elevation model (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     crsapp(S,DEM);  
%
% References
%
%     Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: 
%     uncertainty assessment and smoothing using quantile regression 
%     techniques. Earth Surface Dynamics, 5, 821-839. 
%     [DOI: 10.5194/esurf-5-821-2017]
%
% See also: STREAMobj/crs, STREAMobj/quantcarve, STREAMobj/smooth
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument.')
end

params = struct();

fhandle = figure('Name','crsapp');

%% Controls
hp  = uipanel('FontSize',12,...
             'BackgroundColor','white',...
             'Position',[0 0 0.3 1]);
mt  = uicontrol(hp,'Style','text',...
             'String','controls',...
             'HorizontalAlignment','left',...
             'String','Controls',...
             'Units','normalized',...
             'FontWeight','bold',...
             'Position',[0 .96 1 .04]);
%% K         
s = sprintf('Larger K will increasingly smooth the profiles.');      
cKtext = uicontrol(hp,'Style','text',...
                'HorizontalAlignment','left',...
                'String','K = 1',...
                'Units','normalized',...
                'Position',[0 .85 1 .05],...
                'TooltipString',s);
cK = uicontrol(hp,'Style','slider',...
                'Units','normalized',...
                'Position',[0 .9 1 .045],...
                'Min',-1,...
                'Max',4,...
                'Value',0,...
                'SliderStep',[0.01 0.1],...
                'Callback',@display_Kslider_value,...
                'TooltipString',s);       

%% Gradient            
s = sprintf('Activate downstream monotonicity constraint.');              
cng = uicontrol(hp,'Style','checkbox',...
                'String','Activate downstream monotonicity constraint',...
                'Units','normalized',...
                'Position',[0 .8 1 .045],...
                'TooltipString',s,...
                'Value',1,...
                'Callback',@activateGradientMonotinicity);               
s = sprintf('Imposes a minimum downstream gradient to the profile.');               
cmgtext = uicontrol(hp,'Style','text',...
                'HorizontalAlignment','left',...
                'String','Minimum gradient = 0',...
                'Units','normalized',...                
                'Position',[0 .7 1 .05],...
                'TooltipString',s);
cmg = uicontrol(hp,'Style','slider',...
                'Units','normalized',...
                'Min',0,...
                'Max',0.04,...
                'Value',0,...
                'SliderStep',[0.0001 0.001],...
                'Position',[0 .75 1 .045],...
                'Callback',@display_MGslider_value,...
                'TooltipString',s);              
%% Tau         
s = sprintf('Choose Tau (Quantile between 0 and 1).');      
cTautext = uicontrol(hp,'Style','text',...
                'HorizontalAlignment','left',...
                'String','Tau = 0.5',...
                'Units','normalized',...
                'Position',[0 .6 1 .05],...
                'TooltipString',s);
cTau = uicontrol(hp,'Style','slider',...
                'Units','normalized',...
                'Position',[0 .65 1 .045],...
                'Min',0.0001,...
                'Max',0.9999,...
                'Value',0.5,...
                'SliderStep',[0.001 0.1],...
                'Callback',@display_Tauslider_value,...
                'TooltipString',s);            

s = sprintf('Remove curvature penalty at tributary junctions.');              
cnst = uicontrol(hp,'Style','checkbox',...
                'String','Knicks at tributary junctions',...
                'Units','normalized',...
                'Position',[0 .55 1 .045],...
                'TooltipString',s);            
         
%% Controls                       
% Push button       
s = sprintf('Push to smooth profile.');  
pbh = uicontrol(hp,'Style','pushbutton','String','Calculate now',...
                'Units','normalized',...
                'Position',[0 .40 1 .05],...
                'Callback',@smoothnetwork,...
                'TooltipString',s);               
%% Plot
ha = uipanel('FontSize',12,...
             'BackgroundColor','white',...
             'Position',[.3 0 0.7 1]);
ax = axes('Parent',ha);
box(ax,'on');
hold(ax,'on')
hlorig = plotdz(S,z,'color',[.4 .4 .4]);

hlsmooth = plot(0,min(z));
zs = [];

%% Push
% Push button
s = sprintf('Export smoothed elevations to the base workspace.');  
pexp = uicontrol(hp,'Style','pushbutton','String','Export node attribute list',...
                'Units','normalized',...
                'Position',[0 .35 1 .05],...
                'Callback',@exportnal,...
                'Enable','off',...
                'TooltipString',s);  
% Push button            
s = sprintf('Export parameter settings to the base workspace.'); 
pexpp = uicontrol(hp,'Style','pushbutton','String','Export parameters',...
                'Units','normalized',...
                'Position',[0 .30 1 .05],...
                'Callback',@exportparams,...
                'Enable','off',...
                'TooltipString',s);  

function display_MGslider_value(hObject,callbackdata)
   newval = num2str(hObject.Value);
   set(cmgtext,'String',['Minimum gradient = ' newval]);
end
function display_Kslider_value(hObject,callbackdata)
   newval = num2str(10.^(hObject.Value));
   set(cKtext,'String',['K = ' newval]);
end
function display_Tauslider_value(hObject,callbackdata)
   newval = num2str(hObject.Value);
   set(cTautext,'String',['Tau = ' newval]);
end
function activateGradientMonotinicity(hObject,callbackdata)
    newval = hObject.Value;
    if ~newval
        set(cmg,'Enable','off');
    else
        set(cmg,'Enable','on');
    end
end
function smoothnetwork(hObject,callbackdata)
    set(pbh,'Enable','off','String','Please wait ...');
    set(pexp,'Enable','off')
    set(pexpp,'Enable','off')
    drawnow
    
    params.nonstifftribs = get(cnst,'Value');
    params.K             = 10.^get(cK,'Value');
    params.mingradient   = get(cmg,'Value');
    params.Tau           = get(cTau,'Value');
    
    if ~get(cng,'Value')
        params.mingradient = nan;
    end
    
    assignin('base','params',params);
    
    zs = crs(S,z,params,'split',isempty(gcp('nocreate')));
    if exist('hlsmooth','var')
        delete(hlsmooth)
    end
    hlsmooth = plotdz(S,zs,'color',[1 0 0]);
    drawnow
    
    set(pbh,'Enable','on','String','Calculate now');
    set(pexp,'Enable','on')
    set(pexpp,'Enable','on')
end

function exportnal(hObject,callbackdata)
        
        if ~isempty(zs)
            
            prompt = {'Enter variable name:'};
            title = 'Export';
            lines = 1;
            def = {'zs'};
            answer = inputdlg(prompt, title, lines, def);
            if ~isempty(answer) && isvarname(answer{1})
                assignin('base',answer{1},zs);
            else
                return
            end
        else
            warndlg('No node attribute list available for export.');
        end
end

function exportparams(hObject,callbackdata)
        if isfield(params,'Tau')
            
            
            prompt = {'Enter variable name:'};
            title = 'Export';
            lines = 1;
            def = {'p'};
            answer = inputdlg(prompt, title, lines, def);
            if ~isempty(answer) && isvarname(answer{1})
                assignin('base',answer{1},params);
            else
                return
            end
        else
            warndlg('No node parameters available for export.');
        end
end
end