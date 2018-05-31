function app = manageobjects(hObject,eventdata,app)
% MANAGEOBJECTS provides several simple tools to manage, plot, or export
% topoapp objects.
%
% See also: listfigure, getobjectitem, STEAMobj/modify, 
% STREAMobj2mapstruct, shapewrite (Mapping toolbox), streamorder,
% GRIDobj/crop
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    listicon = imread('listicon.png','png');
    
    % Set-up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',listicon,'TooltipString','Manage objects',...
        'ClickedCallback',{@manageobjects,app});
    
else % execute tool
    
    % Initialize basic figure
    [hf,posYlast] = listfigure(app,app.objects.classes);
    % Add tool-specific options
    addoptions(hf,posYlast,app);
    
end

    %------------------------------------------------------------
    function addoptions(hfig,posYlast,app)

        figure(hfig);
        %--------------------------------------------------------
        % positioning:
        dLs = 0.025; dLl = 0.2833;
        Lft1 = dLs; Lft2 = Lft1+(dLl+dLs+dLs); Lft3 = Lft2+(dLl+dLs+dLs);
        dLftedit1 = 0.12;  dLftedit2 = 0.2133;
        Hpb = 0.06; Hedit = 0.045; Htxt = 0.045;
        Wpb = dLl; Wedit = 0.07; Wtxts = 0.1; Wtxtl = 0.17;
        %--------------------------------------------------------

        %% Export to Google Earth
        %--------------------------------------------------------
        uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.875,Wpb,Hpb],...
            'String','Export to Google Earth',...
            'Callback',{@googleearthexport,app});
        %--------------------------------------------------------

        %% Export to Shapefile
        %--------------------------------------------------------
        uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.8,Wpb,Hpb],...
            'String','Export to Shapefile',...
            'Callback',{@shapefileexport,app});
        %--------------------------------------------------------
        
        %% Plot object
        %--------------------------------------------------------
        userdata1.classes = {'STREAMobj'};
        uih(1) = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.725,Wpb,Hpb],...
            'String','Plot STREAMobj',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
        userdata2.classes = {'WATERSHEDobj'};
        uih(2) = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.725,Wpb,Hpb],...
            'String','Plot WATERSHEDobj',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
        userdata3.classes = {'REACHobj'};
        uih(3) = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.725,Wpb,Hpb],...
            'String','Plot REACHobj',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
        userdata4.classes = {'PROFILEobj'};
        uih(4) = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.725,Wpb,Hpb],...
            'String','Plot PROFILEobj',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
        userdata5.classes = {'SWATHobj'};
        uih(5) = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.725,Wpb,Hpb],...
            'String','Plot SWATHobj',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
            % parameter controls
            
            posY = posYlast;
            h(1) = uicontrol('Style','checkbox',...
                'units','normalized','position',[Lft3,posYlast,Wpb,Htxt],...
                'String','Merge axes',...
                'TooltipString','Plots all items in one figure. Only effective for profile-view');
                posY = posY - 0.05;
            h(2) = uicontrol('Style','checkbox',...
                'units','normalized','position',[Lft3,posY,Wpb,Htxt],...
                'String','Link distance to STREAMobj',...
                'TooltipString','Only effective in profile-view',...
                'Callback',{@link2streamobj,app});
                posY = posY - 0.05; 
            h(3) = uicontrol('Style','edit',...
                'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
                'String','min');
            h(4) = uicontrol('Style','edit',...
                'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
                'String','max');
            h(5) = uicontrol('Style','text',...
                'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
                'String','Strahler','Userdata',h([3,4]));
                posY = posY - 0.05;
            h(6) = uicontrol('Style','checkbox',...
                'units','normalized','position',[Lft3,posY,Wpb,Htxt],...
                'String','Highlight trunkstream');
            h(7) = uicontrol('Style','checkbox',...
                'units','normalized','position',[Lft3,posY,Wpb,Htxt],...
                'Value',1,...
                'String','Left');
                posY = posY - 0.05;
            h(8) = uicontrol('Style','checkbox',...
                'units','normalized','position',[Lft3,posY,Wpb,Htxt],...
                'Value',1,...
                'String','Right');
            % tool options
            h(9) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.275,Wpb,Hpb],...
                'String','Plot profile-view',...
                'Callback',{@plotprofileview,app});
            h(10) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.2,Wpb,Hpb],...
                'String','Plot map-view',...
                'Callback',{@plotmapview,app});
        set(h,'Tag','optionalparameter','Enable','off','Visible','off')
        set(h([1,2,5:8]),'Backgroundcolor',[.8 .8 .8])
        userdata1.handles = h([1:6,9,10]);
        set(uih(1),'Userdata',userdata1,'Visible','off'); % streamobj
        userdata2.handles = h([1,9,10]);
        set(uih(2),'Userdata',userdata2,'Visible','off'); % watersheds
        userdata3.handles = h([1,2,9,10]);
        set(uih(3),'Userdata',userdata3,'Visible','off'); % reaches
        userdata4.handles = h([1,9,10]);
        set(uih(4),'Userdata',userdata4,'Visible','off'); % profiles
        userdata5.handles = h([1,2,7:10]);
        set(uih(5),'Userdata',userdata5,'Visible','off'); clear h; % swathobj
            
        %--------------------------------------------------------
        
        %% Catchment statistics
        %--------------------------------------------------------
        userdata.classes = {'WATERSHEDobj'};
        uih = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.65,Wpb,Hpb],...
            'String','Catchment statistics',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
            % parameter controls...none
            % tool options
            h(1) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.275,Wpb,Hpb],...
                'String','Hypsometry',...
                'Callback',{@catchmentstat,app});
            h(2) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.2,Wpb,Hpb],...
                'String','Basic',...
                'Callback',{@catchmentstat,app});
        userdata.handles = h;
        set(h,'Tag','optionalparameter','Enable','off','Visible','off')
        set(uih,'Userdata',userdata,'Visible','off'); clear h;
        %--------------------------------------------------------

        %% Crop DEM
        %--------------------------------------------------------
        userdata.classes = {'WATERSHEDobj'};
        uih = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.575,Wpb,Hpb],...
            'String','Crop DEM',...
            'Tag','optionaltool',...
            'Callback',{@updatecallbacks});
            % parameter controls
            posY = posYlast; ct = 1;
            h(ct) = uicontrol('Style','text',...
                'units','normalized','position',[Lft3,posY,.15,Htxt],...
                'String','Fill value',...
                'TooltipString','Leave blank for original grid value'); ct=ct+1;
            h(ct) = uicontrol('Style','edit',...
                'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
                'String','NaN',...
                'TooltipString','Leave blank for original grid value'); ct=ct+1;
            % tool options
            h(ct) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.275,Wpb,Hpb],...
                'String','Plot DEM',...
                'Userdata',findobj('Style','edit','TooltipString','Leave blank for original grid value'),...
                'Callback',{@cropdem,app}); ct=ct+1;
            h(ct) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.2,Wpb,Hpb],...
                'String','Export GeoTIFF',...
                'Userdata',findobj('Style','edit','TooltipString','Leave blank for original grid value'),...
                'Callback',{@cropdem,app}); ct=ct+1;
            h(ct) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.125,Wpb,Hpb],...
                'String','Plot surface',...
                'Callback',{@cropdem,app});
        userdata.handles = h;
        set(h,'Tag','optionalparameter','Enable','off','Visible','off')
        set(h(1:2),'Backgroundcolor',[.8 .8 .8])
        set(uih,'Userdata',userdata,'Visible','off'); clear h;
        %--------------------------------------------------------

        %% Export outlets
        %--------------------------------------------------------
        userdata.classes = {'WATERSHEDobj'};
        userdata.handles = [];
        uih = uicontrol('Style','togglebutton',...
            'units','normalized','position',[Lft2,.5,Wpb,Hpb],...
            'String','Export outlets lat,lon',...
            'Tag','optionaltool',...
            'Userdata',userdata,...
            'Callback',{@updatecallbacks});
            % tool options
            ct = 1;
            h(ct) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.275,Wpb,Hpb],...
                'String','Lon/Lat',...
                'Userdata',findobj('Style','edit','TooltipString','Leave blank for original grid value'),...
                'Callback',{@outletexport,app}); ct=ct+1;
            h(ct) = uicontrol('Style','pushbutton',...
                'units','normalized','position',[Lft2,.2,Wpb,Hpb],...
                'String','Map units',...
                'Userdata',findobj('Style','edit','TooltipString','Leave blank for original grid value'),...
                'Callback',{@outletexport,app});
        userdata.handles = h;
        set(h,'Tag','optionalparameter','Enable','off','Visible','off')
        set(uih,'Userdata',userdata,'Visible','off'); clear h;
        %--------------------------------------------------------


    end  % 
    %------------------------------------------------------------
    function updatecallbacks(hobject,~)
        % Display optional tools

        % deactivate optional tools
        h = findobj(gcf,'Style','togglebutton','Value',1);
        h = setdiff(h,hobject);
        if ~isempty(h)
            set(h,'Value',0);
        end

        % remove active parameter controls
        h = findobj(gcf,'Tag','optionalparameter');
        if ~isempty(h)
            set(h,'Visible','off','Enable','off');
        end
        
        userdata = get(hobject,'Userdata');
        if get(hobject,'Value')
            set(userdata.handles,'Visible','on','Enable','on')
        else
            set(userdata.handles,'Visible','off','Enable','off')
        end
    end % 
    %------------------------------------------------------------
    
    %------------------------------------------------------------
    function link2streamobj(~,~,app)
        if get(gcbo,'Value')
            if ~isempty(app.objects.STREAMobj.data)
                % select streamobj
                strobjname = app.objects.STREAMobj.names;
                set(0,'DefaultUIcontrolUnits','pixel');
                [s,ok] = listdlg('PromptString','Select STREAM object to link with:',...
                                'SelectionMode','single',...
                                'ListString',strobjname);
                set(0,'DefaultUIcontrolUnits','normalized');
                if ok
                    set(gcbo,'Userdata',s);
                else
                    set(gcbo,'Value',0)
                end
            else
                fprintf(1,'No STREAMobj to link with.\n');
                set(gcbo,'Value',0)
            end
        else
            set(gcbo,'Userdata',[]);
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
    function [offset] = adjustdistance(classname,ix,app)
        offset = 0;            
        h = findobj('String','Link distance to STREAMobj','Visible','on');
        if get(h,'Value') & ( strcmp(classname,'STREAMobj') || ...
            strcmp(classname,'REACHobj') )
        
            iy = get(h,'Userdata');
            d = app.objects.STREAMobj.data{iy}.distance;
            A = app.objects.STREAMobj.data{iy}.IXgrid(d==0);
            linkobjname = app.objects.STREAMobj.names{iy};
            if ~isempty(A)
                % link distance to lowest point of streamobj
                itemname = app.objects.(classname).names{ix};
                x = app.objects.(classname).data{ix}.x;
                y = app.objects.(classname).data{ix}.y;
                iz = coord2ind(app.DEM,x,y);
                iz = iz(~isnan(iz));
                iz_fa = app.FA.Z(iz);
                iz = iz(iz_fa==max(iz_fa));
                [B,dist] = flowpathextract(app.FD,iz(1));
                [~,~,ib] = intersect(A,B,'stable');
                if ~isempty(ib)
                    offset = dist(ib);
                    fprintf(1,'Distance of %s linked to lowest point of %s.\n',itemname,linkobjname);
                else
                    fprintf(1,'Could not link %s to STREAMobj %s.\n',itemname,linkobjname);
                end
            end
        end
    end % 
    %------------------------------------------------------------
    
    %------------------------------------------------------------
    function googleearthexport(~,~,app)
        if ~exist('kml','file')
            error('Need the kml toolbox by Fernandes de Oliveira (downloadable from Matlab Central)')
        else
            tic
            fprintf(1,'GoogleEarth export in progress...');
            [objectname,~,ix] = getobjectitem(app);
            if ~isempty(ix)
                switch objectname 
                    case 'STREAMobj'
                        k = kml('topoapp_streamobj_kml');
                        for i = 1 : length(ix)
                            [GS,~,~] = STREAMobj2mapstruct(app.objects.(objectname).data{ix(i)});
                            for j = 1 : length(GS)
                                name = [app.objects.(objectname).names{ix(i)},' - ',num2str(j)];
                                [lat,lon] = minvtran(app.DEM.georef.mstruct,GS(j).X,GS(j).Y);
                                k.plot(lon(~isnan(lon)),lat(~isnan(lat)),'name',name);
                            end
                        end
                    case 'SWATHobj'
                        k = kml('topoapp_swathobj_kml');
                        for i = 1 : length(ix)
                            UTMX = app.objects.(objectname).data{ix(i)}.X;
                            UTMY = app.objects.(objectname).data{ix(i)}.Y;
                            for j = 1 : length(UTMX)
                                utmx = reshape(UTMX{j},numel(UTMX{j}),1);
                                utmy = reshape(UTMY{j},numel(UTMY{j}),1);
                                [lat,lon] = minvtran(app.DEM.georef.mstruct,utmx,utmy);
                                k.scatter(lon,lat,'iconColor',[1 0 0 1],'iconScale',0.2);
                            end
                        end
                        
                    otherwise
                        k = kml('topoapp_kml');
                        for i = 1 : length(ix)
                            name = app.objects.(objectname).names{ix(i)};
                            utmx = app.objects.(objectname).data{ix(i)}.x;
                            utmy = app.objects.(objectname).data{ix(i)}.y;
                            [lat,lon] = minvtran(app.DEM.georef.mstruct,utmx',utmy');
                            k.plot(lon,lat,'name',name);
                        end
                end
                fprintf(1,'done. '); toc
                k.run;
            end
        end
    end % 
    %------------------------------------------------------------
    function shapefileexport(~,~,app)
        
        [objectname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            switch objectname 
                case 'STREAMobj'
                    C = [];
                    for i = 1 : length(ix)
                        [tC,x,y] = STREAMobj2mapstruct(app.objects.(objectname).data{ix(i)});
                        jx = coord2ind(app.DEM,x,y);
                        inan = [0;find(isnan(jx))];
                        for k = 1 : length(tC);
                            p1 = inan(k)+1;
                            p2 = inan(k+1)-1;
                            thisjx = jx(p1:p2);
                            maxfa = max(app.FA.Z(thisjx))*app.DEM.cellsize*app.DEM.cellsize;
                            tC(k).MaxArea = maxfa;
                        end
                        C = [C tC]; 
                    end

                otherwise
                    ct = 1;
                    for i = 1 : length(ix)
                        x = app.objects.(objectname).data{ix(i)}.x;
                        y = app.objects.(objectname).data{ix(i)}.y;
%                                 iz = sub2ind(app.DEM.size,y,x);
                        C(ct).Geometry = 'Line';
                        C(ct).BoundingBox = [min(x),min(y); max(x), max(y)];
                        C(ct).X = x;
                        C(ct).Y = y;
                        C(ct).Name = app.objects.(objectname).names{ix(i)};
%                                 C(ct).MaxArea = max(app.FA.Z(iz))*app.DEM.cellsize*app.DEM.cellsize;
                        ct = ct+1;
                    end
            end
            if strcmp('WATERSHEDobj',objectname)
                ct = 1;
                for i = 1 : length(ix)
                    IX = app.objects.(objectname).data{ix(i)}.ix;
                    z = double(app.DEM.Z(IX));
                    g = double(app.G.Z(IX));
                    C(ct).Geometry = 'Polygon';
                    C(ct).Area = length(IX)*app.DEM.cellsize*app.DEM.cellsize;
                    C(ct).MinZ = min(z);
                    C(ct).MaxZ = max(z);
                    C(ct).MeanZ = mean(z);
                    C(ct).MeanSlope = mean(g);
                    ct = ct+1;
                end
            end
            
            if exist('shapewrite','file') % save using mapping toolbox function
                answer = inputdlg('Enter filename (without extension)','Export shapefile');
                if ~isempty(answer)
                    shapefilename = [answer{1},'.shp'];
                    fprintf(1,'*S* Writing file %s...',shapefilename);
                    shapewrite(C,shapefilename);
                    fprintf(1,'done.\n');
                end
                
            else % save as ASCII
                set(0,'DefaultUIcontrolUnits','pixels');
                choice = questdlg('Mapping Toolbox not found. Create ASCII file instead?', ...
                    'Export to Shapefile','Yes','No','Yes');
                set(0,'DefaultUIcontrolUnits','normalized');
                if choice
                    [FileName,PathName,~] = uiputfile('*.txt');
                    fid = fopen([PathName,FileName],'w');
                    
                    fields = setdiff(fieldnames(C),{'BoundingBox'});
                    fprintf(1,'%s,',fields{:});
                    fields = setdiff(fieldnames(C),{'X','Y','BoundingBox'});
                    for i = 1 : length(C)
                        % still work to do
                        
                    end
                    fclose(fid);
                    
                else
                    return;
                end
                
            end
            
        else
            warning('Found no selected object. No data exported.')
        end

    end % 
    %------------------------------------------------------------
    function plotprofileview(~,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            mergeax = get(findobj('String','Merge axes','Visible','on'),'Value');
            if mergeax, figure; end
            for i = 1 : length(ix)
                if ~mergeax, figure; else hold on; end
                distadjust = adjustdistance(classname,ix(i),app);
                switch classname
                    case 'PROFILEobj'
                        plotzobj(app,app.objects.(classname).data{ix(i)},true,distadjust);
                    case 'REACHobj'
                        plotzobj(app,app.objects.(classname).data{ix(i)},false,distadjust);
                    case 'WATERSHEDobj'
                        plotzobj(app,app.objects.(classname).data{ix(i)},false,0);
                    case 'STREAMobj'
                        SObj = app.objects.(classname).data{ix(i)};
                        [SObj] = limit2strahler(SObj);
                        if ~isempty(SObj.ix)
                            hold on
                            plotdz(SObj,app.DEM,'Color',[.7 .7 .7],'doffset',distadjust);
                            plotdz(SObj,app.DEMc,'Color','k','doffset',distadjust);%,'smooth',false);
                            %plotdz(SObj,app.DEMc,'Color','r','doffset',distadjust,'smooth',true,'kernelsize',11);
                            if get(findobj('String','Highlight trunkstream'),'Value')
                                ht = plotdz(trunk(SObj),app.DEMc,'Color','r','doffset',distadjust,'smooth',true,'kernelsize',11);
                                set(ht,'LineWidth',1.5)
                            end
                            hold off
                        end
                    case 'SWATHobj'
                        h = findobj('String','Left','Visible','on');
                        if get(h,'Value'); left=true; else left=false; end
                        h = findobj('String','Right','Visible','on');
                        if get(h,'Value'), right=true; else right=false; end
                        plotdz(app.objects.(classname).data{ix(i)},'left',left,'right',right,'distadjust',distadjust);
                end
                if ~mergeax, updateaxis(gca); title(app.objects.(classname).names{ix(i)}); end
            end
            updateaxis(gca)
        else error('No object to plot')
        end

    end % 
    %------------------------------------------------------------
    function plotmapview(~,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            mergeax = get(findobj('String','Merge axes','Visible','on'),'Value');
            if mergeax, figure; end
            for i = 1 : length(ix)
                if ~mergeax, figure; else hold on; end
                switch classname
                    case 'STREAMobj'
                        SObj = app.objects.(classname).data{ix(i)};
                        [SObj] = limit2strahler(SObj);
                        plotstreamorder(SObj); hold on
                    case 'WATERSHEDobj'
                        x = app.objects.(classname).data{ix(i)}.x;
                        y = app.objects.(classname).data{ix(i)}.y;
                        plot(x,y,'k-'); hold on
                        outlet_ix = app.objects.(classname).data{ix(i)}.outlet;
                        [outr,outc] = ind2sub(app.DEM.size,outlet_ix);
                        plot(app.X(outc),app.Y(outr),'ro');
                        text(app.X(outc),app.Y(outr),app.objects.(classname).names{ix(i)});
                    case 'SWATHobj'
                        SWObj = app.objects.(classname).data{ix(i)};
                        h = findobj('String','Left','Visible','on');
                        if get(h,'Value'); left=true; else left=false; end
                        h = findobj('String','Right','Visible','on');
                        if get(h,'Value'), right=true; else right=false; end
                        SWObj.plot('left',left,'right',right);
                    otherwise
                        x = app.objects.(classname).data{ix(i)}.x;
                        y = app.objects.(classname).data{ix(i)}.y;
                        plot(x,y,'k-'); hold on
                        text(x(1),y(1),app.objects.(classname).names{ix(i)});
                end
                if ~mergeax, updateaxis(gca); end
            end
            axis image
            updateaxis(gca)
        else error('No object to plot')
        end
    end % 
    %------------------------------------------------------------
    function catchmentstat(~,~,app)
    % Plot some simple catchment statistics
        [~,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            switch get(gcbo,'String')
                case 'Basic'
                    % make buffer around streams
                    Sgrid = STREAMobj2GRIDobj(app.S);
                    SE = strel('disk',3,4);
                    Sgrid.Z = imdilate(Sgrid.Z,SE);
                    INO = find(Sgrid.Z==1);

                    fprintf(1,'\nName\t\t\tArea(km^2)\tMinZ(m)\tMaxZ(m)\tMeanZ(m)\tMeanSlope(deg)\tMeanSlopeBuff(deg)\n');
                    fprintf(1,'%s\n',repmat('-',1,100));
                    for i = 1:length(ix)
                        IX = app.objects.WATERSHEDobj.data{ix(i)}.ix;
                        z = app.DEM.Z(IX);
                        g = app.G.Z(IX);
                        IY = setdiff(IX,INO);
                        gb = app.G.Z(IY);

                        fprintf(1,'%s\t\t%1.2f\t\t%i\t\t%i\t%i\t\t%1.1f\t%1.1f\n',...
                            app.objects.WATERSHEDobj.names{ix(i)},...
                            numel(z)*app.DEM.cellsize*app.DEM.cellsize/1e6,...
                            min(z),max(z),round(mean(z)),mean(g),mean(gb));
                    end
                    fprintf(1,'\n');

                case 'Hypsometry'
                    figure
                    for i = 1:length(ix)
                        IX = app.objects.WATERSHEDobj.data{ix(i)}.ix;
                        tDEM = crop(app.DEM,IX);
                        [rf,elev] = hypscurve(tDEM);
                        h(i) = plot(rf,elev,'-','Color',rand(1,3));
                        hold on
                    end
                    hold off
                    xlabel('Area (%)')
                    ylabel('Elevation (m)')
                    wsnames = app.objects.WATERSHEDobj.names(ix);
                    legend(h,wsnames)
                    updateaxis(gca)
            end
        end

    end % 
    %------------------------------------------------------------
    function cropdem(hobject,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            IX = [];
            for i = 1 : length(ix)
                IX = [IX;app.objects.(classname).data{ix(i)}.ix]; % indices to watershed pixels
            end
            h = get(hobject,'Userdata');
            fillvalue = str2double(get(h,'String'));
            M = app.DEM;
            M.Z = 0;
            M.Z(IX)=1;
            DEMcropped = crop(app.DEM,M,fillvalue);
            caller = get(hobject,'String');
            switch caller
                case 'Plot DEM'
                    figure, imageschs(DEMcropped);
                    waitfor(h)
                case 'Plot surface'
                    figure, surf(DEMcropped);
                    waitfor(h)
                case 'Export GeoTIFF'
                    GRIDobj2geotiff(DEMcropped)
            end
        end
    end % 
    %------------------------------------------------------------
    function outletexport(~,~,app)
    % Export outlet locations as lat,lon coordinate list
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            IX = nan(size(ix));
            for i = 1 : length(ix)
                IX(i) = app.objects.(classname).data{ix(i)}.outlet;
            end
            [x,y] = ind2coord(app.DEM,IX);
            if strcmp(get(gcbo,'String'),'Lon/Lat')
                fprintf(1,'Watershed\tLon\t\t\tLat\n');
                [y,x] = projinv(app.DEM.georef,x,y);
            else
                fprintf(1,'Watershed\tX\t\t\tY\n');
            end
            wsnames = app.objects.(classname).names(ix);
            for i = 1 : length(x)
            fprintf(1,'%s\t%f\t%f\n',wsnames{i},x(i),y(i));
            end
        else
            fprintf(1,'No watersheds selected.\n');
        end
    end % 
    %------------------------------------------------------------
    
end % 
   
    