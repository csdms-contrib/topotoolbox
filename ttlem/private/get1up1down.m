function [ii,kk] = get1up1down(i,k,A)
% get receiver of receiver and giver of giver...
% ii and kk will have nans where neighbors do not exist, % either because, there is none % or because the upstream node has a smaller % contribution area.

i = double(i);
k = double(k);

nrc = numel(A);

% Downstream neighbor
kk = nan(nrc,1);
kk(i) = k;
kk(i) = kk(k);
kk = kk(i);

% Upstream neighbor
ii = accumarray(k,i,[nrc 1],@getindexwithlargerarea,nan); ii = ii(i);


function ix = getindexwithlargerarea(ix)

[~,mix] = max(A(ix));
ix = ix(mix);
end

end