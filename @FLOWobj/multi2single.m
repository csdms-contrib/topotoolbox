function [FD,S] = multi2single(FD,varargin)
%MULTI2SINGLE converts multiple to single flow direction
%
% Syntax
%
%     FD = multi2single(FDm)
%     FD = multi2single(FDm,pn,pv,...)
%     [FD,S] = multi2single(FDm,pn,pv,...)
%
% Description
%
%     multi2single converts a multiple flow direction FLOWobj FDm into a 
%     single flow direction FLOWobj. Additional arguments enable to convert
%     only the channelized part defined by a minimum upstream area from
%     multi to single. 
%
% Input arguments
%
%     FDm   multiple flow directions
%     
%     Parameter name/value pairs
%
%     'minarea'      minimum upstream area required to initiate streams. The
%                    default is 0. Higher values will result in flow
%                    networks that have multiple flow directions up to the
%                    defined upstream area, and single flow directions
%                    downstream.
%     'unit'         unit of value determined with the parameter 'minarea':
%                    'pixels' (default) or 'mapunits'.
%     'channelheads' linear index into DEM with channelheads. 
%
% Output arguments
%
%     FD    FLOWobj. If 'minarea' == 0 (default), then FD will be a FLOWobj  
%           with single flow directions. Otherwise, FLOWobj will be of type 
%           'multi'.
%     S     STREAMobj (only applicable if 'minarea' is set >0 or
%           channelheads are provided.)
%     
%
% Example
% 
%       DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%       FD  = FLOWobj(DEM,'preprocess','carve','mex',true);
%       DEM = imposemin(FD,DEM,0.0001);
%       S   = STREAMobj(FD,'minarea',1000);
%       DEMc = DEM;
%       DEMc.Z(S.IXgrid) = DEMc.Z(S.IXgrid)-100; 
%       FD  = FLOWobj(DEMc,'multi');
%       FD  = multi2single(FD,'minarea',100);
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. May, 2019

p = inputParser;
p.FunctionName = 'FLOWobj/multi2single';
addParameter(p,'minarea',0,@(x) isscalar(x));
addParameter(p,'unit','pixels',@(x) ischar(validatestring(x, ...
                            {'pixels', 'mapunits'})));
addParameter(p,'channelheads',[])
                        
parse(p,varargin{:});

switch FD.type
    case 'single'
        
        % do nothing
        
    otherwise
        
        if p.Results.minarea > 0 || ~isempty(p.Results.channelheads)
            
            if isempty(p.Results.channelheads)
                
                unit    = validatestring(p.Results.unit,{'pixels', 'mapunits'});
                switch unit
                    case 'mapunits'
                        minarea = p.Results.minarea/(FD.cellsize.^2);
                    otherwise
                        minarea = p.Results.minarea;
                end
                
                % Find stream network initiation points. This should not be
                % done via thresholding flow accumulation derived from multiple
                % flow directions because we may erroneously include pixels in
                % diverging parts of the river network.
                
                ix  = FD.ix;
                ixc = FD.ixc;
                fr  = FD.fraction;
                A   = ones(FD.size);
                ini = false(FD.size);
                for r = 1:numel(ix)
                    if A(ix(r)) < minarea
                        A(ixc(r)) = A(ixc(r)) + fr(r)*A(ix(r));
                    elseif A(ix(r)) >=minarea
                        ini(ix(r)) = true;
                    end
                end
                
            else
                if isa(p.Results.channelheads,'GRIDobj')
                    ini = p.Results.channelheads.Z > 0;
                else
                    ini = false(FD.size);
                    ini(p.Results.channelheads) = true;
                end
            end
            
            if nargout == 2
                channelheads = find(ini);
            end
            
        end
        
        RR = (1:numel(FD.ix))';
        IX = double(FD.ix);
        S  = sparse(RR,IX,FD.fraction,max(RR),max(IX));
        [~,ii] = max(S,[],1);
        I  = false(size(RR));
        I(ii) = true;
        
        if p.Results.minarea <= 0 && isempty(p.Results.channelheads)
        
            FD.ix  = FD.ix(I);
            FD.ixc = FD.ixc(I);
            FD.fraction = [];
        
            FD.type = 'single';
            
            if nargout == 2
                S = [];
            end
            
        else
            
            ix = FD.ix(I);
            ixc = FD.ixc(I);
            
            for r = 1:numel(ix)
                ini(ixc(r)) = ini(ixc(r)) || ini(ix(r));
            end
            
            I = I | ~ini(FD.ix);
            
            FD.fraction(ini(FD.ix)) = 1;
            FD.ix = FD.ix(I);
            FD.ixc = FD.ixc(I);
            FD.fraction = FD.fraction(I);
            
            FD.type = 'multi';
            
            if nargout == 2
                FD.type = 'single';
                S = STREAMobj(FD,'channelheads',channelheads);
                FD.type = 'multi';
            end
        end
            
        
end
end