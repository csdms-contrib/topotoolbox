function [H,Pc] = dbentropy(M,siz,ix)

% entropy of drainage basin delineation
%
% Syntax
%
%     [H,S] = dbentropy(M,siz,ix)
%
% Description
%
%     dbentropy calculates the Shannon Entropy (H) of a digital elevation
%     model. H quantifies the uncertainty associated with each pixel in the
%     DEM to drain towards a specific outlet.
%
% Input parameters
%
%     M     multiple flow direction matrix
%     ix    linear indices of outlets. Outlets must not be located in the
%           upslope flow network of another outlet.
%     siz   size of the digital elevation model
%
% Output arguments
%  
%     H     Shannon Entropy
%     S     structure array with numel(ix) elements with additional 
%           information to each outlet. S.P contains the probability grid
%           for each outlet. S.IX is the outlet's linear index. CA01-99 are
%           the error bounds for the drainage area.
%
% Example
%
%     load exampledem
%     M = flowdir(X,Y,fillsinks(dem));
%     [H,S] = dbentropy(M,size(dem));
%     imageschs(X,Y,dem,H)
%     hold on
%     plot(X([S.IX]),Y([S.IX]),'kx')
%
% See also: FLOWDIR, DRAINAGEBASINS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 27. October, 2011

if nargin == 2;
    ix = find(((sum(M,2)) == 0) & (sum(M,1)' > 0));
    nrix = numel(ix);
    % check the number of outlets and break the problems into chunks if
    % necessary
    
    chunksize = 10;
    
    if nrix > chunksize;
        chunks = true;
    else
        chunks = false;
    end
    
    % In addition, suggest to output only some of the larger fuzzy drainage basins 
    if nargout >= 2 && chunks;
        
        
        mess  = ['The number of probability matrices is quite large (k = ' num2str(nrix) ').\n' ...
                 'It is suggested to return only the probability matrices of\n'...
                 'the outlets with the largest drainage basins. Please provide\n'...
                 'a number of how many you want. The default is n = ' num2str(chunksize) '. n : '];
        pout = input(mess);
        if isempty(pout)
            pout = 10;
        end
        
        A = (speye(prod(siz))-M')\ones(prod(siz),1);
        [~,ixsort] = sort(A(ix),'descend');
        ix = ix(ixsort);
        
    end
    nrc  = prod(siz);
    
else
    nrc  = prod(siz);
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

if nargout == 2;
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
        
        for r = 1:nrix;
            Pc(r).P = full(reshape(P(:,r),siz));
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
            
        end
    end
    
    sP    = sum(P,2);
    sP    = sparse(max(1-sP,0));

    H = reshape(-full(sum(spfun(@(x) x.*log2(x),[P sP]),2)),siz);

    
else
    H = zeros(siz);
    counter = 1;
    sP = zeros(prod(siz),1);
    
    for r = 1:ceil(nrix/chunksize);
        cols = ((r-1)*chunksize + 1):r*chunksize;
        cols(cols>nrix) = [];
        
        P = (speye(nrc)-M)\S(:,cols);
        sP = sum(P,2)+sP;
        H = H + reshape(-full(sum(spfun(@(x) x.*log2(x),P),2)),siz);

        if nargout >= 2
            
            for rr = cols(1):cols(end);
                
                if rr<=pout
                    Pc(rr).P = full(reshape(P(:,counter),siz));
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
                if counter >= chunksize;
                    counter = 1;
                else
                    counter = counter+1;
                end
            
                
            end
        end
    end 
%     H = H + reshape(1 - sP.*log2(sP),siz);
end

H = max(min(H,1),0);


% plotit = true;
% 
% if plotit
%     r = numel(Pc);
%     bar(1:r,[Pc.CA50],...
%         'BaseValue',0,...
%         'EdgeColor','none',...
%         'FaceColor',[.6 .6 .6]);
%     hold on
%     plot(repmat(1:r,2,1),[[Pc.CA01];[Pc.CA99]],'k','LineWidth',2);
%     hold off
%     
% end



% function to calculate the log of an arbitrary base
% function y = logb(x,b)
% 
% if b == 1;
%     y = log(x);
% elseif b == 2;
%     y = log2(x);
% elseif b == 10;
%     y = log10(x);
% else
%     y = log(x)./log(b);
% end

