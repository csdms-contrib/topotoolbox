function [CS,locS,order] = STREAMobj2cell(S,ref,n)

%STREAMOBJ2CELL convert instance of STREAMobj to cell array of stream objects
%
% Syntax
%
%     CS = STREAMobj2cell(S)
%     CS = STREAMobj2cell(S,ref)
%     CS = STREAMobj2cell(S,'outlets',n)
%     [CS,locS] = ...
%     [CS,locS,order] = STREAMobj2cell(S,'tributaries');
%
% Description
%
%     STREAMobj2cell divides a STREAMobj into a number of STREAMobj stored
%     in a cell array. 
%
%     STREAMobj2cell(S,'outlets') divides a STREAMobj into its strongly
%     connected components. This means that individual STREAMobjs are
%     derived as individual trees of the stream network, i.e. individual
%     drainage basins. This is the default. In this case, CS has as many
%     elements as there are outlets in the stream network.
%     
%     STREAMobj2cell(S,'channelheads') derives individual STREAMobj as
%     single streams emanating from each channelhead. In this case, CS has
%     as many elements as there are channelheads in the stream network.
%
%     STREAMobj2cell(S,'tributaries') derives individual STREAMobj as
%     tributaries. A stream is a tributary until it reaches a stream with a 
%     longer downstream flow distance.     
%     
% Input arguments
%
%     S     instance of STREAMobj
%     ref   reference for deriving individual STREAMobj. Either 'outlets'
%           (default) or 'channelheads', or 'tributaries'
%     n     if ref is 'outlets', n determines the number of n largest trees
%           to be placed in elements of CS.
%
% Output arguments
%
%     CS    cell array with elements of CS being instances of STREAMobj
%     locS  cell array with linear indices into node attribute lists of S
%     order output argument only if ref is set to 'tributaries'. Vector
%           with numel(CS) elements where each element refers to an order 
%           of the tributaries. 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     z   = getnal(S,DEM);
%     [CS,locS] = STREAMobj2cell(S);
%     plotdz(CS{21},z(locS{21}))
%
%
% See also: FLOWobj2cell
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. March, 2017

if nargin == 1
    ref = 'outlets';
    getall = true;
    n   = inf;
elseif nargin == 2
    ref = validatestring(ref,{'outlets','channelheads','tributaries'},'STREAMobj2cell','ref',2);
    getall = true;
    n   = inf;
elseif nargin == 3
    ref = validatestring(ref,{'outlets'},'STREAMobj2cell','ref',2);
    validateattributes(n,{'numeric'},{'>',1},'STREAMobj2cell','n',3);
    getall = false;
end

switch lower(ref)
    case 'outlets'
        
        nrc = numel(S.x);
        M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
        
        [~,p,~,r] = dmperm(M | M' | speye(nrc));
               
        nc = length(r) - 1;
        if getall || nc < n
            % label matrix
            L = zeros(nrc,1);
            for tt = 1:nc
                L(p(r(tt):r(tt+1)-1)) = tt;
            end
        else
            nc = n;
            L = zeros(nrc,1)+nc;
            
            [~,dd] = sort(diff(r),'descend');
            
            counter = 1;

            for tt = dd(1:nc) %1:min(nc,k);
                L(p(r(tt):r(tt+1)-1)) = counter;
                counter = counter + 1;
            end
            
        end
        
        % put each individual tree into an own element in CS.
        CS = cell(1,nc);
        
        % adapt new STREAMobj to the reduced network
        LL = L;
        
        if nargout == 2
            locS = cell(1,nc);
        end
        
        for r = 1:nc
            CS{r} = S;
            L     = LL==r;
            I     = L(CS{r}.ix);
            CS{r}.ix  = CS{r}.ix(I);
            CS{r}.ixc = CS{r}.ixc(I);
            
            IX    = cumsum(L);
            CS{r}.ix  = IX(CS{r}.ix);
            CS{r}.ixc = IX(CS{r}.ixc);
            
            CS{r}.x   = CS{r}.x(L);
            CS{r}.y   = CS{r}.y(L);
            CS{r}.IXgrid   = CS{r}.IXgrid(L);
            
            if nargout == 2
                locS{r} = find(L);
            end
                
            
        end
        
        
        
    case 'channelheads'
        ch     = streampoi(S,'channelheads','logical');
        ixcix  = zeros(numel(S.IXgrid),1);
        ixcix(S.ix) = 1:numel(S.ix);
        
        ixchannel = find(ch);
        nc = numel(ixchannel);
        
        CS = cell(1,nc);
        if nargout == 2
            locS = cell(1,nc);
        end
        for r = 1:nc
            IX = ixchannel(r);
        
            c = 1;
            while ixcix(IX) ~= 0
                c = c+1;
                IX(c) = S.ixc(ixcix(IX(end)));
            end
            
            L = false(size(S.IXgrid));
            L(IX) = true;
            
            % adapt new STREAMobj to the reduced network
            
            CS{r} = S;
            I     = L(CS{r}.ix);
            CS{r}.ix  = CS{r}.ix(I);
            CS{r}.ixc = CS{r}.ixc(I);
            
            IX    = cumsum(L);
            CS{r}.ix  = IX(CS{r}.ix);
            CS{r}.ixc = IX(CS{r}.ixc);
            
            CS{r}.x   = CS{r}.x(L);
            CS{r}.y   = CS{r}.y(L);
            CS{r}.IXgrid   = CS{r}.IXgrid(L);
            
            if nargout == 2
                locS{r} = find(L);
            end
            
            
            
        end
        
    case 'tributaries'
        if nargout < 3
            CS = tributaries(S);
        else
            CS = tributariesinclorder(S);
            order = CS(2:2:end);
            order = horzcat(order{:});
            CS = CS(1:2:end);
        end
        
        if nargout >= 2
            locS = cell(size(CS));
            for r = 1:numel(CS)
                [~,locS{r}] = ismember(CS{r}.IXgrid,S.IXgrid);
            end
        end
        
        return
        

        
end

% check for validity of Ss
valid = true(1,nc);
for r = 1:nc
    if numel(CS{r}.x) == 1
        valid(r) = false;
    end
end
CS = CS(valid);

if nargout == 2
    locS = locS(valid);
end
end


%% Recursively scan for tributaries
function Ctribs = tributaries(S)
% Recursive extraction of tributaries


C = STREAMobj2cell(S);
Ctribs = cell(0);

for r = 1:numel(C)
    St = trunk(C{r});
    S2 = modify(C{r},'tributaryto2',St);
    if isempty(S2.ix)
        % do nothing
        Ctribs   = [Ctribs {St}];
    else
        Ctribs = [Ctribs {St} tributaries(S2)];
    end
end
end



function Ctribs = tributariesinclorder(S,orderin)
% Recursive extraction of tributaries with order


C = STREAMobj2cell(S);
Ctribs = cell(0);

if nargin == 1
    orderin = 1;
end
    

for r = 1:numel(C)
    St = trunk(C{r});
    S2 = modify(C{r},'tributaryto2',St);
    if isempty(S2.ix)
        % do nothing
        Ctribs   = [Ctribs {St orderin}];
    else
        Ctribs = [Ctribs {St orderin} tributariesinclorder(S2, orderin+1)];
    end
end
end

