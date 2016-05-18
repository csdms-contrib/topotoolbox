function CS = STREAMobj2cell(S,ref,n)

% convert instance of STREAMobj to cell array of stream objects
%
% Syntax
%
%     CS = STREAMobj2cell(S)
%     CS = STREAMobj2cell(S,ref)
%     CS = STREAMobj2cell(S,ref,n)
%
% Description
%
%     STREAMobj2cell places individual STREAMobj in elements of a cell
%     array. These individual STREAMobjs may either be derived as
%     individual trees of the stream network, e.g. individual drainage
%     basins. This is the default. In this case, CS has as many elements as
%     there are outlets in the stream network. Individual STREAMobjs can
%     also be derived as single streams emanating from each channelhead. In
%     this case, CS has as many elements as there are channelheads in the
%     stream network.
%     
%     
% Input arguments
%
%     S     instance of STREAMobj
%     ref   reference for deriving individual STREAMobj. Either 'outlets'
%           (default) or 'channelheads'
%     n     if ref is 'outlets', n determines the number of n largest trees
%           to be placed in elements of CS.
%
% Output arguments
%
%     CS    cell array with elements of CS being instances of STREAMobj
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 3. May, 2013

if nargin == 1;
    ref = 'outlets';
    getall = true;
    n   = inf;
elseif nargin == 2;
    ref = validatestring(ref,{'outlets','channelheads'},'STREAMobj2cell','ref',2);
    getall = true;
    n   = inf;
elseif nargin == 3;
    ref = validatestring(ref,{'outlets','channelheads'},'STREAMobj2cell','ref',2);
    validateattributes(n,{'numeric'},{'>',1},'STREAMobj2cell','n',3);
    getall = false;
end

switch lower(ref)
    case 'outlets'
        
        nrc = numel(S.x);
        M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
        
        [~,p,~,r] = dmperm(M | M' | speye(nrc));
               
        nc = length(r) - 1;
        if getall || nc < n;
            % label matrix
            L = zeros(nrc,1);
            for tt = 1:nc;
                L(p(r(tt):r(tt+1)-1)) = tt;
            end
        else
            nc = n;
            L = zeros(nrc,1)+nc;
            
            [~,dd] = sort(diff(r),'descend');
            
            counter = 1;

            for tt = dd(1:nc); %1:min(nc,k);
                L(p(r(tt):r(tt+1)-1)) = counter;
                counter = counter + 1;
            end
            
        end
        
        % put each individual tree into an own element in CS.
        CS = cell(1,nc);
        
        % adapt new STREAMobj to the reduced network
        LL = L;
        for r = 1:nc;
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
        end
        
        
        
    case 'channelheads'
        ch     = streampoi(S,'channelheads','logical');
        ixcix  = zeros(numel(S.IXgrid),1);
        ixcix(S.ix) = 1:numel(S.ix);
        
        ixchannel = find(ch);
        nc = numel(ixchannel);
        
        CS = cell(1,nc);
        
        for r = 1:nc;
            IX = ixchannel(r);
        
            c = 1;
            while ixcix(IX) ~= 0
                c = c+1;
                IX(c) = S.ixc(ixcix(IX(end)));
            end
            
            L = false(size(S.IXgrid));
            L(IX) = true;
            
            % adapt new STREAMobj to the reduced network
%             LL = L;
            
                CS{r} = S;
%                 L     = LL==r;
                I     = L(CS{r}.ix);
                CS{r}.ix  = CS{r}.ix(I);
                CS{r}.ixc = CS{r}.ixc(I);

                IX    = cumsum(L);
                CS{r}.ix  = IX(CS{r}.ix);
                CS{r}.ixc = IX(CS{r}.ixc);

                CS{r}.x   = CS{r}.x(L);
                CS{r}.y   = CS{r}.y(L);
                CS{r}.IXgrid   = CS{r}.IXgrid(L);
            
        end
        
end

% check for validity of Ss
valid = true(1,nc);
for r = 1:nc;
    if numel(CS{r}.x) == 1;
        valid(r) = false;
    end
end
CS = CS(valid);

