function [varargout] = sort(DIN)
%SORT   Sort divide segments by network structure.
%
% Syntax
%
%     DOUT = sort(DIN)
%
% Description
%
%     SORT changes the order of divide segments in a divide
%     object, starting from the endpoints and moving up the
%     divide network. Each segment will be ordered such that
%     more external nodes are higher up in the list of linear
%     indices IX.
%
% Input
%
%     DIN       instance of class DIVIDEobj
%
% Output
%
%     DOUT      instance of class DIVIDEobj
%
% Example
%
%     DOUT = sort(DIN);
%
% See also: DIVIDEobj, DIVIDEobj/divnet
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: Nov 2018

%  input properties: D.IX, D.ep, D.jct, D.jctedg
% output properties: D.IX, D.issorted

DOUT = DIN;

% Prepare
ixep = DIN.ep;
ixjct = DIN.jct;
nedges = DIN.jctedg;

M = onl2struct(DIN.IX);
[M.flag] = deal(true);

% Divide classification loop
segments = struct;
jctcnt = nedges-1; % number of segments before junction becomes an edge

% % indicate here endorheic basins?

jcthit = zeros(size(ixjct)); % count times a jct was met during loop

ct = 0;
%fprintf(1,'Divide classification:\n');
while ~isempty(ixep)
    
    %ct
    % find segments with active endpoints
    ia = ismember(vertcat(M.st),ixep);
    [locr,locc] = find(ia);
    
    % flip segments
    ix = find(locc==2);
    for i = 1 : length(ix)
        thisr = locr(ix(i));
        M(thisr).IX = [flipud(M(thisr).IX(1:end-1));NaN];
        M(thisr).st = fliplr(M(thisr).st);
    end
    
    % loop over active endpoints
    for i = 1 : length(locr)
        thisr = locr(i);
        if M(thisr).flag
            newep = M(thisr).st(2); % second, because flipped
            addjct = ismember(ixjct,newep);
            if nnz(addjct)>0 % label junction
                jcthit(addjct) = jcthit(addjct)+1;
            end
            % transfer segment
            ct = ct+1;
            segments(ct).ix = M(thisr).IX;
            M(thisr).flag = false;
        end
    end % for-endpoints loop
    
    % find new endpoints
    ix = jcthit>=jctcnt;
    if sum(ix)>0
        % set new endpoints
        ixep = ixjct(ix);
        % remove old endpoints from junctions
        [ixjct,ia] = setdiff(ixjct,ixep,'stable');
        jcthit = jcthit(ia);
        jctcnt = jctcnt(ia);
    else
        ixep = [];
    end
    
end % while-endpoints loop


DOUT.IX = vertcat(segments.ix);
DOUT.issorted = true;

if numel(DIN.IX)~=numel(DOUT.IX)
    warning('Warning: Some divide segments could not be ordered. Results may be erroneous.')
end

varargout{1} = DOUT;
if nargout==2
    varargout{2} = setdiff(DIN.IX,DOUT.IX,'stable');
end


end


