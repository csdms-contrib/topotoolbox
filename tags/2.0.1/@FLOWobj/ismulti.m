function tf = ismulti(FD)
% check if FD is multi or single flow direction

tf = false;
r = 2;
while ~tf && r<=numel(FD.ix);
    tf = FD.ix(r) == FD.ix(r-1);
    r  = r+1;
end
end