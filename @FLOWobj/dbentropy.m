function [OUT,Pc] = dbentropy(FD,ix)

%DBENTROPY entropy of drainage basin delineation
%
% Syntax
%
%     H = dbentropy(FD)
%     [H,P] = dbentropy(FD)
%     [H,P] = dbentropy(FD,IX)
%     [H,P] = dbentropy(FD,S)
%
% Description
%
%     dbentropy calculates the drainage basin entropy H (Schwanghart and
%     Heckmann 2012) based on multiple flow directions stored in FLOWobj
%     FD. Drainage basin entropy is derived for each pixel and quantifies
%     the uncertainty to assign an outlet pixel. Pixel that drain towards a
%     single outlet pixel have low entropy values, whereas pixels with high
%     values could be assigned to two or more outlet pixels. Drainage basin
%     entropy relies on a probabilistic interpretation of multiple flow 
%     directions.
%
% Input parameters
%
%     FD    multiple flow direction FLOWobj
%     IX    linear indices of outlets
%     S     STREAMobj from which linear indices of outlets are derived
%
% Output arguments
%  
%     H     Shannon Entropy grid (GRIDobj)
%     P     structure array with numel(ix) elements with additional 
%           information to each outlet. P.P contains the probability grid
%           for each outlet. P.IX is the outlet's linear index. CA01-99 are
%           the error bounds for the drainage area.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM = fillsinks(DEM);
%     FD = FLOWobj(DEM,'multi');
%     IX = [107807 769639];
%     [H,SB] = dbentropy(FD,IX);
%     imageschs(DEM,H,'colormap',flipud(hot))
%
% Reference
%
%     Schwanghart, W., Heckmann, T. (2012): Fuzzy delineation of drainage 
%     basins through probabilistic interpretation of diverging flow algorithms. 
%     Environmental Modelling and Software, 33, 106-113. 
%     [DOI: 10.1016/j.envsoft.2012.01.016]
%
% See also: drainagebasins, FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. March, 2019

M = FLOWobj2M(FD);
siz = FD.size;

if nargin == 1
    ix = find(((sum(M,2)) == 0) & (sum(M,1)' > 0));
    nrix = numel(ix);
    % check the number of outlets and break the problems into chunks if
    % necessary
    
    chunksize = 100;
    
    if nrix > chunksize
        chunks = true;
    else
        chunks = false;
    end
    
    % In addition, suggest to output only some of the larger fuzzy drainage basins 
    if nargout >= 2 && chunks
              
        mess  = ['The number of probability matrices is quite large (k = ' num2str(nrix) ').\n' ...
                 'It is suggested to return only the probability matrices of\n'...
                 'the outlets with the largest drainage basins. Please provide\n'...
                 'a number of how many you want. The default is n = ' num2str(chunksize) '. n : '];
             
        pout = input(mess);
        
        if isempty(pout)
            pout = 10;
        end
        
        % calculate flow accumulation
        A = flowacc(FD);
        [~,ixsort] = sort(A.Z(ix),'descend');
        ix = ix(ixsort);
        
    end
    nrc  = prod(siz);
    
else
    nrc  = prod(siz);
    if isa(ix,'STREAMobj')
        ix = streampoi(ix,'outlets','ix');
    end
    nrix = numel(ix);
    
    % create sinks at provided indices
    SINKS = ones(nrc,1);
    SINKS(nrix) = 0;
    M    = spdiags(SINKS,0,nrc,nrc)*M;
    
    chunks = false;
end



% S is a matrix where each outlet has its own column. 
S = sparse(ix,1:nrix,1,nrc,nrix);

% B is the membership (partition) matrix (?) (fuzzy partition matrix). 
% B contains the membership grades. (calculated by membership function below)

if nargout >= 2
    Pc = struct('P',{},'IX',{},'CA50',{},'CA99',{},...
                    'CA95',{},'CA66',{},'CA33',{},...
                    'CA05',{},'CA01',{});
    % linear index
    pl = (1:prod(siz))';
    % percentiles
    percentiles = [.99 .95 .66 .50 .33 .05 .01];
end
    

if ~chunks
    % solve system of equations to get probabilities
    P = (speye(nrc)-M)\S;
    
    if nargout >= 2
        
        for r = 1:nrix
            Pc(r).P = GRIDobj(FD);
            Pc(r).P.Z = full(reshape(P(:,r),siz));
            Pc(r).IX = ix(r);
            ptemp = sort(full(P(:,r)),'descend');
            [~,m] = unique(ptemp);
            perc  = interp1(ptemp(m),pl(m),percentiles);
            
            Pc(r).CA50 = perc(4);
            Pc(r).CA99 = perc(1);
            Pc(r).CA95 = perc(2);
            Pc(r).CA66 = perc(3);
            Pc(r).CA33 = perc(5);
            Pc(r).CA05 = perc(6);
            Pc(r).CA01 = perc(7);
            
            Pc(r).percentiles = percentiles;
            Pc(r).area = perc;
            
            
        end
    end
    
    sP    = sum(P,2);
    sP    = sparse(max(1-sP,0));

    H = reshape(-full(sum(spfun(@(x) x.*log2(x),[P sP]),2)),siz);

    
else
    H = zeros(siz);
    counter = 1;
    sP = zeros(prod(siz),1);
	
	h = waitbar(0,'Starting');
	iters = ceil(nrix/chunksize);
    
    for r = 1:iters
		waitbar(r/iters,h,[num2str(r) ' / ' num2str(iters)]);
        cols = ((r-1)*chunksize + 1):r*chunksize;
        cols(cols>nrix) = [];
        
        P = (speye(nrc)-M)\S(:,cols);
        sP = sum(P,2)+sP;
        H = H + reshape(-full(sum(spfun(@(x) x.*log2(x),P),2)),siz);

        if nargout >= 2
            
            for rr = cols(1):cols(end)
                
                if rr<=pout
                    Pc(rr).P = GRIDobj(FD);
                    Pc(rr).P.Z = full(reshape(P(:,counter),siz));
                    ptemp    = sort(Pc(rr).P(:),'descend');
                else
                    Pc(rr).P = [];
                    ptemp    = sort(full(P(:,counter)),'descend');
                end
                
                Pc(rr).IX = ix(rr);
                [~,m] = unique(ptemp);
                perc  = interp1(ptemp(m),pl(m),percentiles);
            
                Pc(rr).CA50 = perc(4);
                Pc(rr).CA99 = perc(1);
                Pc(rr).CA95 = perc(2);
                Pc(rr).CA66 = perc(3);
                Pc(rr).CA33 = perc(5);
                Pc(rr).CA05 = perc(6);
                Pc(rr).CA01 = perc(7);
                
                
                % reset counter if larger than chunksize
                if counter >= chunksize
                    counter = 1;
                else
                    counter = counter+1;
                end
            
                
            end
        end
    end 
	close(h);
%     H = H + reshape(1 - sP.*log2(sP),siz);
end

H = max(min(H,1),0);
OUT = GRIDobj(FD);
OUT.Z = H;
OUT.name = 'fuzzy drainage basins';


