function FD = multi2single(FD)
%MULTI2SINGLE converts multiple to single flow direction
%
% Syntax
%
%     FD = multi2multi(FD)
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016

switch FD.type
    case 'single';
        
        % do nothing
        
    otherwise
        RR = (1:numel(FD.ix))';
        IX = double(FD.ix);
        S  = sparse(RR,IX,FD.fraction,max(RR),max(IX));
        [~,ii] = max(S,[],1);
        I  = false(size(RR));
        I(ii) = true;
        
        FD.ix = FD.ix(I);
        FD.ixc = FD.ixc(I);
        FD.fraction = [];
        
        FD.type = 'single';
        
end
end