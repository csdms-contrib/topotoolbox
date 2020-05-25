function DOUT = divdist(DIN)
% DIVDIST   Assign distance to divide segments
%
% Syntax
%
%     D2 = divdist(D)
%
% Description
%
%     DIVDIST assigns to each segment in the divide network a vector with
%     the maximum distance from an endpoint measured aloing the sorted
%     divide network. The distance is stored in the divide object in the
%     field 'distance'
%
% Input
%
%     D       instance of class DIVIDEobj
%
% Output
%
%     D2      instance of class DIVIDEobj
%
% Example
%
%     D = divdist(D);
%
% See also: DIVIDEobj, DIVIDEobj/sort
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


DOUT = DIN;

% Prepare
M = onl2struct(DIN.IX);
[M.dist] = deal([]);

% note where the end of the segment is the start of a new one
st = vertcat(M.st);
st1 = st(:,1);
st2 = st(:,2);
%[~,lib] = ismember(st2,st1);

% set distance (0) for endpoints
ixep = ismember(st1,DIN.ep);
[M(ixep).dist] = deal(0);
[M.checked] = deal(false);
% segments with two endpoints
is2ep = sum(ismember(st,DIN.ep),2)==2;

for i = 1 : length(M)
    %ixc = lib(i);
    ix = M(i).IX;
    M(i).dist = [max(M(i).dist) + (0:numel(ix)-2)'.*DIN.cellsize; NaN];
    if is2ep(i) % two endpoints
        d = M(i).dist;
        M(i).dist = min([d d([end-1:-1:1 end])],[],2);
    end
    M(i).checked = true;
    
    lib = ismember(st1,st2(i));
    if nnz(lib)>0 %ixc>0
        allixc = find(lib);
        for k = 1 : length(allixc)
            ixc = allixc(k);
            if not(M(ixc).checked)
                M(ixc).dist = [M(ixc).dist; max(M(i).dist)];
            end
        end
    end
end


DOUT.distance = vertcat(M.dist);

end
