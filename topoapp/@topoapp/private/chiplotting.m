function app = chiplotting(hObject,eventdata,app)
% CHIPLOTTING provides access to TopoToolbox function 'chiplot'
%
% See also: topoapp/listfigure, CHIPLOT
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    chiploticon = imread('chiicon.png','png');
    
    % Set-up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',chiploticon,'TooltipString','Chi-plot',...
        'ClickedCallback',{@chiplotting,app});
    
else % execute tool
    
    if isempty(app.S)
        warning('No STREAMobj found. Use FLOW routing button first')
    else
        % Initialize basic figure
        [hf,posYlast] = listfigure(app,{'STREAMobj'});
        % Add tool-specific options
        addoptions(hf,posYlast,app);
    end
end
    %------------------------------------------------------------
    function hfig = addoptions(hfig,posYlast,app)
        
        figure(hfig)
        %--------------------------------------------------------
        % positioning:
        dLs = 0.025; dLl = 0.2833;
        Lft1 = dLs; Lft2 = Lft1+(dLl+dLs+dLs); Lft3 = Lft2+(dLl+dLs+dLs);
        dLftedit1 = 0.12;  dLftedit2 = 0.2133;
        Hpb = 0.06; Hedit = 0.045; Htxt = 0.045;
        Wpb = dLl; Wedit = 0.07; Wtxts = 0.1; Wtxtl = 0.17;
        %--------------------------------------------------------

        %% Parameter controls
        posY = posYlast;
        h(1) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','m/n',...
            'TooltipString','Leave blank for automatic least-squares search');
        h(2) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','',...
            'TooltipString','Leave blank for automatic least-squares search');
        posY = posY - 0.05;
        h(3) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','A0');
        h(4) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','1e6',...
            'TooltipString','Reference area in square map units');
        posY = posY - 0.05;
        h(5) = uicontrol('Style','checkbox',...
            'String','Highlight trunk stream',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',1);
        posY = posY - 0.05;
        h(6) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
            'String','min');
        h(7) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','max');
        h(8) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
            'String','Strahler','Userdata',h([6,7]));

        %% Tool options
        h(9) = uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.275,Wpb,Hpb],...
            'String','Calculate m/n from reach',...
            'Callback',{@calcmnreach,app});
        h(10) = uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.2,Wpb,Hpb],...
            'String','Chi-plot',...
            'Userdata',h([2,4,5]),...
            'Callback',{@mkchiplot,app});
        set(h([1,3,5,8]),'Backgroundcolor',[.8 .8 .8])
        %--------------------------------------------------------
        
    end %
    %------------------------------------------------------------
    function mkchiplot(~,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            for i = 1 : length(ix) % loop over objects
                s = app.objects.(classname).data{ix(i)};
                name = app.objects.(classname).names{ix(i)};
                [s] = limit2strahler(s);
                if isempty(s)
                    fprintf(1,'Current limits imposed by strahler order create empty STREAMobj./n');
                else
                    h = findobj('Style','edit','TooltipString','Leave blank for automatic least-squares search');
                    mn = str2double(get(h,'String'));
                    if isnan(mn); mn = []; end
                    h = findobj('Style','edit','TooltipString','Reference area in square map units');
                    A0 = str2double(get(h,'String'));
                    if isnan(A0); A0 = 1e6; end
                    h = findobj('Style','checkbox','String','Highlight trunk stream');
                    trunkstream = get(h,'Value');
                    if trunkstream
                        C = chiplot(s,app.DEMc,app.FA,'mn',mn,'A0',A0,'plot',true,'mnplot',isempty(mn),'trunkstream',trunk(s));
                    else
                        C = chiplot(s,app.DEMc,app.FA,'mn',mn,'A0',A0,'plot',true,'mnplot',isempty(mn));
                    end
                    updateaxis(gca)
                    title(name)
                    fprintf(1,'%s: ks(m/n=%f) = %f +/- %f\n',name,C.mn,C.ks,C.betase*(C.a0^C.mn));
                end
            end
            h = findobj('Style','edit','TooltipString','Leave blank for automatic least-squares search');
            set(h,'String',num2str(C.mn));
        end
        
    end % 
    %------------------------------------------------------------
    function calcmnreach(~,~,app)
        if ~isempty(app.objects.REACHobj.data)
            set(0,'DefaultUIcontrolUnits','pixel');
            itemnames = app.objects.REACHobj.names;
            [s,ok] = listdlg('PromptString','Select a reach:',...
                'SelectionMode','single',...
                'ListString',itemnames);
            set(0,'DefaultUIcontrolUnits','normalized');
            
            if ok
                ix = app.objects.REACHobj.data{s}.ix;
                W = app.DEM;
                W.Z(:) = false;
                W.Z(ix) = true;
                s = STREAMobj(app.FD,W);
                c = chiplot(s,app.DEM,app.FA,'plot',false,'mnplot',false);
                h = findobj('Style','edit','TooltipString','Leave blank for automatic least-squares search');
                set(h,'String',num2str(c.mn));
            end
            
        else
            fprintf(1,'Found no reaches.\n');
        end
        
    end % 
    %------------------------------------------------------------
    function [SObj] = limit2strahler(SObj)
        hsto = findobj('String','Strahler','Visible','on');
        if ~isempty(hsto)
            hst = get(hsto,'Userdata');
            stmin = str2double(get(hst(1),'String'));
            if isnan(stmin);
                set(hst(1),'String','min');
            else
                SObj = modify(SObj,'streamorder',['>=',num2str(stmin)]);
            end
            stmax = str2double(get(hst(2),'String'));
            if isnan(stmax);
                set(hst(2),'String','max');
            else
                SObj = modify(SObj,'streamorder',['<=',num2str(stmax)]);
            end
        end
    end % 
    %------------------------------------------------------------


end



