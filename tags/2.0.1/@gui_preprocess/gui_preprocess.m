classdef gui_preprocess < handle


properties
    
    % Data properties and display information
    DEM  % digital elevation model (GRIDobj)
    DEMf % filled digital elevation model (GRIDobj)
    X    % coordinate vector
    Y    % "
    DEMundo % digital elevation model for undo operation (GRIDobj)
    RGB
    
    
    % Layout properties
    fh   % figure handle
    menu % menu handle
    tbhb % handle zoom toolbar
    ax   % handle axis for plotting
    im   % handle image
    zm   % handle zoom object
    lh   % vector of handles to line objects
    mpanel % panel for buttons etc.
    helptext % help text
    
    methodgroup
    togglecarve
    togglefill
    
    panelcarve
    carvegroup
    carvemeth
    carvevaltext
    carvevaledit
    bcarve
    
    panelfill
    fillgroup
    fillmeth
    fillvaltext
    fillvaledit
    bfill
    
end

properties (Dependent = true, SetAccess = private)
    SINKDEPTH 
    SINKS % logical array that is true for sinks
    
end



methods
    function this_app = gui_preprocess(DEM)
        % run app
        
        % check input
        
        this_app.DEM  = DEM;
        this_app.DEMf = fillsinks(DEM);
        [this_app.X, this_app.Y] =  getcoordinates(DEM);
        this_app.X = this_app.X(:)';
        this_app.Y = this_app.Y(:);
        
        layout(this_app)
        this_app.RGB = imageschs(this_app.DEM,this_app.SINKDEPTH./(this_app.SINKS));
               
    end
    
    %% Functions for calculating dependent properties
    function SINKDEPTH = get.SINKDEPTH(app)
        SINKDEPTH = app.DEMf-app.DEM;
    end
    function SINKS = get.SINKS(app)
        SINKS = app.DEMf > app.DEM;
    end
    %% On-Close function
    function closefcn(~,~,app)
        delete(app.fh);
        delete(app);
    end

end

end














