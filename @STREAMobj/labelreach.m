function [label,varargout] = labelreach(S,varargin)

%LABELREACH create node-attribute list with labelled reaches
%
% Syntax
%
%     label = labelreach(S)
%     label = labelreach(S,pn,pv,...)
%     [label,ix] = ...
%     [label,x,y] = ...
%
% Description
%
%     labelreach creates a node attribute list where each element refers to
%     a label of a reach. A reach is defined as a section of the river
%     network between two confluences. An additional criteria for labelling
%     can be reach distance so that reaches do not exceed a specified length
%     (see option seglength).
%
% Input arguments
%
%     S       STREAMobj
%     Parameter name value pairs
%     'seglength'   segment (label) length
%     'shuffle'     randomize labels
%
% Output arguments
%
%     label   node-attribute list with labels
%     ix      linear index into GRIDobj of split locations
%     x,y     coordinates of split locations
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     [l,x,y] = labelreach(S,'seglength',1000,'shuffle',true);
%     plotc(S,l)
%     hold on
%     plot(x,y,'k+')
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

p = inputParser;
p.FunctionName = 'STREAMobj/labelreach';
addParamValue(p,'seglength',inf,@(x) isscalar(x) && x>0);
addParamValue(p,'shuffle',false,@(x) isscalar(x));
addParamValue(p,'exact',false,@(x) isscalar(x));
parse(p,varargin{:});


if ~p.Results.exact
    % The first approach will label all nodes in the stream network.
    % It will split reaches according to segment length but also at
    % confluences. This will result in segment lengths that are not exactly
    % the provided value seglength.
    label = streampoi(S,{'bconfl','out'},'logical');
    label = cumsum(label).*label;
    
    ix = S.ix;
    ixc = S.ixc;
    for r = numel(ix):-1:1
        if label(ix(r)) == 0
            label(ix(r)) = label(ixc(r));
        end
    end
    
    % there may be nodes that are neither receivers or givers. These nodes
    % still have the label zero.
    maxlabel = max(label);
    iszero   = label==0;
    label(iszero) = (1:nnz(iszero))+maxlabel;
    
    
    if ~isinf(p.Results.seglength)
        maxdist   = p.Results.seglength;
        cs        = S.cellsize;
        dist      = S.distance;
        label2    = accumarray(label,(1:numel(S.x))',[max(label) 1],@(x) {split(x)});
        labeltemp = zeros(size(label));
        for r = 1:max(label)
            labeltemp(label2{r}(:,2)) = label2{r}(:,1);
        end
        [~,~,label] = unique([label labeltemp],'rows');
    end
else
    % If we choose the exact approach, then labelreach will not split at
    % confluences. Parts close to channelheads that have length less than
    % the value of seglength will get the label zero.
    [CS,locS] = STREAMobj2cell(S,'tributaries');
    labeltrib = zeros(size(S.x));
    for r = 1:numel(locS)
        labeltrib(locS{r}) = r;
    end
    
    label     = zeros(size(S.x));
    Clabel = cell(size(CS));
    
    for iter = 1:numel(CS)
        d = CS{iter}.distance;
        d = mod(d,p.Results.seglength);
        I = double(gradient(CS{iter},d) < 0);
        I(streampoi(CS{iter},'outlet','logical')) = 1;
        
        ix = CS{iter}.ix;
        ixc = CS{iter}.ixc;
        
        for r = numel(ix):-1:1
            I(ix(r)) = I(ixc(r)) + I(ix(r));
        end
        Clabel{iter} = I;
    end
    
    for r = numel(CS):-1:1
        label(locS{r}) = Clabel{r};
    end
    
    % get unique labels
    [~,~,label] = unique([label labeltrib],'rows');
    
end
    


%% Shuffle labels
if p.Results.shuffle
    [uniqueL,~,ix] = unique(label);
    uniqueLS = randperm(numel(uniqueL));
    label = uniqueLS(ix);
end

label = label(:);

if nargout >= 2
    I = label(S.ix)~=label(S.ixc);
    
    if nargout == 2
        varargout{1} = S.IXgrid(S.ixc(I));
    else
        varargout{1} = S.x(S.ixc(I));
        varargout{2} = S.y(S.ixc(I));
    end
end
    

%% Split function
function sp = split(ix)

ix = ix(:);
d = dist(ix)-min(dist(ix));
maxd = max(d);
if maxd < maxdist
    sp = [ones(size(ix)) ix];
else
    nrsegs   = ceil(maxd./maxdist);
    edges    = (0:maxd/nrsegs:(maxd+cs/100));
    [~,~,sp] = histcounts(d,edges);
    sp = [sp(:) ix];
end
end

end 
    
