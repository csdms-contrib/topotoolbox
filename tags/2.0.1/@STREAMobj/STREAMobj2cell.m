function CS = STREAMobj2cell(S)

% convert instance of STREAMobj to cell array of stream objects
%
% Syntax
%
%     CS = STREAMobj2cell(S)
%
% Description
%
%     STREAMobj2cell transfers each individual tree of the stream network 
%     into a separate instance of STREAMobj in a cell array of STREAMobjs.
%
% Input arguments
%
%     S     instance of STREAMobj
%
% Output arguments
%
%     CS    cell array with elements of CS being instances of STREAMobj
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 3. May, 2013

nrc = numel(S.x);
M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);

[~,p,~,r] = dmperm(M | M' | speye(nrc));
nc = length(r) - 1;

% label matrix
L = zeros(nrc,1);
for tt = 1:nc;
    L(p(r(tt):r(tt+1)-1)) = tt;
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


