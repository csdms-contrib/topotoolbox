function carvemanually(varargin)
% Manual carving function
app = varargin{1};

keepgoing = true;
counter = 0;
xx = zeros(2,1);
yy = zeros(2,1);

while keepgoing && get(app.bcarve,'Value')
    k = waitforbuttonpress;
    
    if k == 0
        p = get(app.ax,'CurrentPoint');
        x = p(1);
        y = p(1,2);
        
        if get(app.fh,'CurrentObject') ~= app.im;           
            pause(.1);
            if get(app.bcarve,'Value') == 0;
                 break
            end        
        else
            counter = counter + 1;
            hold on
            app.lh(end+1) = plot(app.ax,x,y,'*k'); %#ok<*AGROW>
            hold off
            xx(counter) = x;
            yy(counter) = y;
            
            if counter == 2;
                hold on
                app.lh(end+1) = plot(xx,yy,'k--');
                hold off
                counter = 0;
                
                IX = coord2ind(app.DEM,xx,yy);
     
                switch find(cell2mat(get(app.carvemeth,'Value')));
                    % linear ramp, minima imposition
                    case {1 2}
                        straightcarving(app,IX);
                    case 3
                        lcpcarving(app,IX);
                        
                end
                
                % write back in matrix
                updateimage(app)
            end
        end
    end
    

end
deletelines(app)


end