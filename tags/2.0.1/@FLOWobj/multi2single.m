function FD = multi2single(FD)
switch FD.type
    case 'single';
        
        % do nothing
        
    otherwise
        nre = numel(FD.ix);
        I   = false(nre,1);
        
        IX = FD.ix(1);
        LargestFraction = FD.fraction(1);
        IXLargestFraction = 1;
        
        
        for r = 2:nre;
            if FD.ix(r) == IX;
                if FD.fraction(r)>=LargestFraction;
                    IXLargestFraction = r;
                end
            else
                I(IXLargestFraction) = true;
                IX = FD.ix(r);
                LargestFraction = FD.fraction(r);
                IXLargestFraction = r;
            end
        end
        I(IXLargestFraction) = true;
        
        
        
        FD.ix = FD.ix(I);
        FD.ixc = FD.ixc(I);
        FD.fraction = [];
        
        FD.type = 'single';
        
end
end