function fillmanually(app)
% Manual filling function
keepgoing = true;

% label matrix
L = bwlabel(app.SINKS.Z);

while keepgoing && get(app.bfill,'Value')
    
    try
        k = waitforbuttonpress;
        
        if k == 0
            p = get(app.ax,'currentpoint');
            x = p(1);
            y = p(1,2);
            if get(app.fh,'CurrentObject') ~= app.im;
                
                pause(.1);
                if get(app.bfill,'Value') == 0;
                    break
                end
            else
                IX = coord2ind(app.X,app.Y,x,y);
                I  = L==L(IX);
                minz = min(min(app.DEM.Z(imdilate(I,ones(3)) & ~I)));
                app.DEM.Z(I) = minz;
                updateimage(app);
                
            end
        end
    catch
        return
    end
end
end