function A = flowacc_lm(dem,W,varargin)

% flow accumulation (upslope area) for LARGE Digital Elevation Models
%
% Syntax
%
%     A = flowacc_lm(dem)
%     A = flowacc_lm(dem,W0)
%     A = flowacc_lm(dem,[],propertyname,propertyvalue,...)
%     A = flowacc_lm(dem,W0,propertyname,propertyvalue,...)
%
% Description
%
%     Multiple flowdirection and flowaccumulation algorithm that routes 
%     through flat terrain. The function is designed to handle LARGE 
%     matrices. This function is memory efficient by looping through tiles
%     of the digital elevation model (DEM). It does not return the flow
%     direction matrix, since it requires relatively large amounts of
%     memory, in particular when the multiple flow direction algorithm is
%     used.
%
% Input
%
%     dem       digital elevation model
%     W0        weight matrix same size as dem 
%               (default: ones(size(dem)))
%
% Properties
%
%     propertyname     propertyvalues/description
%
%     'type'          'multi' (default): multiple flowdirection (dinf)
%                     'single': single flow direction (d8). Flow 
%                     occurs only along the steepest descent
%
%     'mode'          'default': deterministic flow
%                     'random': totally random flow to downward 
%                     neighbors
%                     'randomized': deterministic flow with noise
%
%     'exponent'       exponent governing the relation between flow
%                      direction and slope. Default is 1, which means,  
%                      there is a linear relation. You may want to 
%                      increase the exponent when flow direction should 
%                      rather follow a steepest descent (single) flow 
%                      direction (e.g. 5). This option is only effective 
%                      for multiple flowdirection.
%
%     'fillsinks'     'no' (default): pits are not removed
%                     'yes': pits are filled 
%
%     'routeflats'    'single' (default), 'multi' or 'none': determines
%                     the way, the algorithm routes through flats
%
%     'waitmessage'   'on' (default) displays information on the progress 
%                     of calculation in the command window. 'off' does 
%                     not.
%
%     'tilesize'      50000 (default). 
%
%
% Output
%
%     A         flowaccumulation (upslope area) grid 
%
% Example
% 
%     load exampleDEM
%     % upscale DEM by linear interpolation
%     dem = interp2(dem,3);
%     X   = interp2(X,3);
%     Y   = interp2(Y,3);
%     A   = flowacc_lm(dem);
% 
% 
% 
%
% Acknowledgements:
% flowacc_lm uses three subfunctions. Two of them are available on
% the file exchange and have been copied here. 
% 1.  parseargs by Malcolm Wood 
%     http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?object
%     Id=10670&objectType=file
% 2.  wdisp by us
%     http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?object
%     Id=1436&objectType=file
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 8. July, 2009



ini_W = 1;

if nargin == 1;
    % initial amount of water in each cell
    W = ones(size(dem)+2)*ini_W;
else
    if isempty(W);
        W = ones(size(dem)+2)*ini_W; 
    else
        % check size 
        if ~isequal(size(W),size(dem));
            error('TopoToolbox:incorrectinput',...
                  'W0 must have same size as the DEM')
        end
        W = padarray(W,[1 1],0);
    end
end

% check input using PARSEARGS
params.type        = {'multi','single'};
params.mode        = {'default','randomized','random'};
params.routeflats  = {'route','cross','none'};
params.waitmessage = {'on','off'};
params.fillsinks   = {'no','yes'};
params.tilesize    = 50000;
params.exponent    = 1;
params = parseargs(params,varargin{:});


% wait message
switch params.waitmessage
    case 'on'
        fprintf('  process:\n')
end

% pad array
padval     = min(dem(:)-1);
dem        = padarray(dem,[1 1],padval);

% check if there are any nans
inan       = isnan(dem);
if any(inan(:));
    flagnan = true;
    dem(inan) = padval;
else
    flagnan = false;
    clear inan
end

% fillsinks
switch params.fillsinks
    case 'yes'
        % fill sinks
        dem = fillsinks(dem);
        switch params.waitmessage
            case 'on'
                fprintf('     - sinks filled \n')
        end
end

% route through flats
switch params.routeflats
    case 'route'
        [icf,icdf] = routeflats(dem,params.type);
        switch params.waitmessage
            case 'on'
                fprintf('     - routed through flats \n')
        end
    case 'cross'
        [icf,icdf,W0] = crossflats(dem,params.type);
        switch params.waitmessage
            case 'on'
                fprintf('     - crossed flats \n')
        end
        W(icf) = W0;
    otherwise
        icf = [];
        icdf = [];      
end


% iteratively calculate upslope area
I = identifyflats(dem);
I = I(:);

% sort dem
[IXs,IXs] = sort(dem(:),'descend'); 
if flagnan
    IXs(inan(IXs)) = [];
end

% remove indices on the edge of the dem;
Ipad = false(size(dem));
Ipad(2:end-1,2:end-1) = true;
IXs = IXs(Ipad(IXs));

IXs_isflat = I(IXs);

% return chunk lengths (should more or less have the length of
% params.tilesize. Yet, this is not always possible since plateaus
% at one specific elevation have to be taken all at one step.

% general values
siz      = size(dem);
nIXs     = numel(IXs);
nrtiles  = ceil(nIXs/params.tilesize);

% initialize area
A   = W;
ix1 = 1;
ix2 = params.tilesize;
maxtilesize = 360000; % take care with adjusting this value


% waitmessage
switch params.waitmessage
    case 'on'
        fprintf('     - calculate upslope area\n');
end

% waitbar
% percent per run

for run = 1:nrtiles;
    % adjust ix2 if located in a plateau
    ix2 = min(ix2,nIXs);
        
    
    if IXs_isflat(ix2) && ~ix2==nIXs;
        % find next non plateau cell
        ix2t = find(~IXs_isflat(ix2:end),1,'first')+ix2-1;
        % may this new index exceed maxtilesize
        if ix2t-ix1+1 > maxtilesize;
            % if yes look in the different direction
            ix2t = find(~IXs_isflat(ix1:ix2),1,'last')+ix1-1;
            
            if ix2t == ix1;
                error('TopoToolbox:plateau_problem',...
                  'Plateau sizes are too large too handle')
            end
        end
        ix2 = ix2t;
    end
    
    % waitmessage
    switch params.waitmessage
        case 'on'
            p = ix2/nIXs*100;
            tl=wdisp( 0,sprintf('       percent done %.2f',p));
    end
    
    IX = IXs(ix1:ix2);
    If = IXs_isflat(ix1:ix2);
    
    % calculate slopes in all downslope directions
    [ic,icd,slope] = calcslope(IX(~If),dem);
    III = ~isnan(slope);
    ic  = ic(III);
    icd = icd(III);
    slope   = slope(III);
    
    III = ismember(icf,IX(If));
    ic  = [[ic;icf(III)] [icd;icdf(III)]];
    s   = [slope;ones(sum(III),1)];
    
    
    [IXu,m,n]  = unique(ic);
    n          = reshape(n,[],2);
    
    At         = A(IXu);    
    nrAt       = numel(At);
    
    M          = sparse(n(:,1),n(:,2),s,nrAt,nrAt);

    % ******************************************************************
    % Randomization of amount transferred to another cell
    switch params.mode;
        case 'random'
            M = abs(sprandn(M));
        case 'randomized'
            % randomize coefficient. The higher, the more random
            rc = 0.01;
            M = M + (rc * abs(sprandn(M)));
        otherwise
    end
    % ******************************************************************
    % single flow direction
    switch params.type
        case 'single'
            [m,IX2] = max(M,[],2);
            i = m==0;
            IX1 = (1:nrAt)';
            IX1(i) = [];
            IX2(i) = [];
            M = sparse(IX1,IX2,ones(length(IX1),1),nrAt,nrAt);
        otherwise
            if params.exponent ~= 1;
                M = M.^params.exponent;
            end
    end

    % ******************************************************************
    % Row standardization of M only necessary when multiple flow dir
    switch params.type
        case 'multi'
            s     = full(sum(M,2));
            ir    = s ~= 0;
            s(ir) = 1./s(ir); 
            M     = spdiags(s,0,nrAt,nrAt) * M;
        otherwise
    end
    
    % calculate drainage area
    At = (speye(nrAt) - M')\At;
    
    A(IXu) = At;
    
    switch params.waitmessage
        case 'on'
            tl=wdisp(tl); %#ok
    end
    
    
    ix1 = ix2+1;
    ix2 = ix1+params.tilesize;

end
    
% reshape to orginal size    
A = reshape(A,siz);

% depad array
A = A(2:end-1,2:end-1);

% and this is it...   









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IXflats = unique(icf);
%        
% % flats have been removed now
% 
% 
% % flow routing from here on
% % sort indices of dem according to their elevation from up to down 
% [IXs,IXs] = sort(dem(:),'descend'); 
% % remove indices on the padded rim around the dem, remove indices of flats
% % and nans
% I4                  = true(size(dem));
% I4(2:end-1,2:end-1) = false;
% I4(IXflats)         = false;
% if flagnan
%     I4(inan)  = false;
% end
% IXs(I4(IXs))  = [];
% clear I4
% 
% % indices that have been indicated as flats are also removed
% 
% % nr of runs
% nIXs     = numel(IXs);
% 
% siz      = size(dem);
% nrc      = prod(siz);
% 
% % initial upslope area
% if isempty(W);
%     A = ones(nrc,1);
% else
%     W = padarray(W,[1 1],0);
%     W = W(:);
%     if numel(W) == numel(dem);
%         A = W;
%     else
%         error('The weight matrix must have same size as the DEM')
%     end
% end
% 
% % waitbar
% % percent per run
% nrrun    = ceil(nIXs/params.tilesize);
% nrlast   = mod(nIXs,params.tilesize);
% 
% switch params.waitmessage
%     case 'on'
%         fprintf('     - calculate upslope area\n');
% end
% 
% % loop through tiles
% for run = 1:nrrun;
%     % waitbar
%     switch params.waitmessage
%         case 'on'
%             p = run/nrrun*100;
%             tl=wdisp( 0,sprintf('       percent done %.2f',p));
%     end
%     if run<nrrun || nrlast==0;
%         IX = IXs((run-1)*params.tilesize + 1 : run*params.tilesize);
%     else
%         IX = IXs((run-1)*params.tilesize + 1 : nIXs);
%     end
%         
%     % calculate slopes in all downslope directions
%     [ic,icd,slope] = calcslope(IX,dem);
%     
%     if flagnan
%         I = inan(icd);
%         ic(I)    = [];
%         icd(I)   = [];
%         slope(I) = [];
%     end
%     
%     % enqueue flats
%     I          = dem(icf) <= max(dem(IX)) & dem(icf) >= min(dem(IX));
%     
%     ic         = [[ic;icf(I)] [icd;icdf(I)]];
%     s          = [slope;ones(sum(I),1)];
%     
%     [IXu,m,n]  = unique(ic);
%     n          = reshape(n,[],2);
%     
%     At         = A(IXu);    
%     nrAt       = numel(At);
%     
%     M          = sparse(n(:,1),n(:,2),s,nrAt,nrAt);
%     
%     icf(I) = [];
%     icdf(I) = [];
% 
%     % ******************************************************************
%     % Randomization of amount transferred to another cell
%     switch params.mode;
%         case 'random'
%             M = abs(sprandn(M));
%         case 'randomized'
%             % randomize coefficient. The higher, the more random
%             rc = 0.01;
%             M = M + (rc * abs(sprandn(M)));
%         otherwise
%     end
%     % ******************************************************************
%     % single flow direction
%     switch params.type
%         case 'single'
%             [m,IX2] = max(M,[],2);
%             i = m==0;
%             IX1 = (1:nrAt)';
%             IX1(i) = [];
%             IX2(i) = [];
%             M = sparse(IX1,IX2,ones(length(IX1),1),nrAt,nrAt);
%         otherwise
%     end
% 
%     % ******************************************************************
%     % Row standardization of M only necessary when multiple flow dir
%     switch params.type
%         case 'multi'
%             s     = full(sum(M,2));
%             ir    = s ~= 0;
%             s(ir) = 1./s(ir); 
%             M     = spdiags(s,0,nrAt,nrAt) * M;
%         otherwise
%     end
%     
%     % calculate drainage area
%     At = (speye(nrAt) - M')\At;
%     
%     A(IXu) = At;
%     
%     switch params.waitmessage
%         case 'on'
%             tl=wdisp(tl);
%     end
% 
% end
% 
% A = reshape(A,siz);
% 
% % depad array
% A(:,1)   = [];
% A(:,end) = [];
% A(1,:)   = [];
% A(end,:) = [];
% 
% % and this is it...




% subfuntion slope calculation

function [ic,icd,s] = calcslope(IX,dem)

siz = size(dem);
n   = numel(IX);
icd = zeros(n,8,'uint32');
s   = zeros(n,8);
cs  = 1;
csd = hypot(cs,cs);

% right
icd(:,1) = IX + siz(1);
s(:,1)   = dem(IX)-dem(icd(:,1));
% right up
icd(:,2) = icd(:,1)+1;
s(:,2)   = (dem(IX)-dem(icd(:,2)))./csd;
% right down
icd(:,3) = icd(:,1)-1;
s(:,3)   = (dem(IX)-dem(icd(:,3)))./csd;

% left
icd(:,4) = IX - siz(1);
s(:,4)   = dem(IX)-dem(icd(:,4));
% left up
icd(:,5) = icd(:,4)+1;
s(:,5)   = (dem(IX)-dem(icd(:,5)))./csd;
% right down
icd(:,6) = icd(:,4)-1;
s(:,6)   = (dem(IX)-dem(icd(:,6)))./csd;

% up
icd(:,7) = IX - 1;
s(:,7)   = dem(IX)-dem(icd(:,7));
% down
icd(:,8) = IX + 1;
s(:,8)   = dem(IX)-dem(icd(:,8));

ic  = repmat(IX,8,1);
icd = icd(:);
s   = s(:);

I     = s<0;
s(I)  = [];
ic(I) = [];
icd(I)= [];









% subfunction parseargs

function X = parseargs(X,varargin)

%PARSEARGS - Parses name-value pairs
%
% Behaves like setfield, but accepts multiple name-value pairs and provides
% some additional features:
% 1) If any field of X is an cell-array of strings, it can only be set to
%    one of those strings.  If no value is specified for that field, the
%    first string is selected.
% 2) Where the field is not empty, its data type cannot be changed
% 3) Where the field contains a scalar, its size cannot be changed.
%
% X = parseargs(X,name1,value1,name2,value2,...) 
%
% Intended for use as an argument parser for functions which multiple options.
% Example usage:
%
% function my_function(varargin)
%   X.StartValue = 0;
%   X.StopOnError = false;
%   X.SolverType = {'fixedstep','variablestep'};
%   X.OutputFile = 'out.txt';
%   X = parseargs(X,varargin{:});
%
% Then call (e.g.):
%
% my_function('OutputFile','out2.txt','SolverType','variablestep');

% The various #ok comments below are to stop MLint complaining about
% inefficient usage.  In all cases, the inefficient usage (of error, getfield, 
% setfield and find) is used to ensure compatibility with earlier versions
% of MATLAB.

remaining = nargin-1; % number of arguments other than X
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
% Take input arguments two at a time until we run out.
while remaining>=2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname,fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X,fieldname); %#ok
        newvalue = varargin{count+1};
        if iscell(oldvalue)
            % Cell arrays must contain strings, and the new value must be
            % a string which appears in the list.
            if ~iscellstr(oldvalue)
                error(sprintf('All allowed values for "%s" must be strings',fieldname));  %#ok
            end
            if ~ischar(newvalue)
                error(sprintf('New value for "%s" must be a string',fieldname));  %#ok
            end
            if isempty(find(strcmp(oldvalue,newvalue))) %#ok
                error(sprintf('"%s" is not allowed for field "%s"',newvalue,fieldname));  %#ok
            end
        elseif ~isempty(oldvalue)
            % The caller isn't allowed to change the data type of a non-empty property,
            % and scalars must remain as scalars.
            if ~strcmp(class(oldvalue),class(newvalue))
                error(sprintf('Cannot change class of field "%s" from "%s" to "%s"',...
                    fieldname,class(oldvalue),class(newvalue))); %#ok
            elseif numel(oldvalue)==1 & numel(newvalue)~=1 %#ok
                error(sprintf('New value for "%s" must be a scalar',fieldname));  %#ok
            end
        end
        X = setfield(X,fieldname,newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
% Check that we had a value for every name.
if remaining~=0
    error('Odd number of arguments supplied.  Name-value pairs required');
end

% Now find cell arrays which were not modified by the above process, and select
% the first string.
notmodified = find(~modified);
for i=1:length(notmodified)
    fieldname = fields{notmodified(i)};
    oldvalue = getfield(X,fieldname); %#ok
    if iscell(oldvalue)
        if ~iscellstr(oldvalue)
            error(sprintf('All allowed values for "%s" must be strings',fieldname)); %#ok
        elseif isempty(oldvalue)
            error(sprintf('Empty cell array not allowed for field "%s"',fieldname)); %#ok
        end
        X = setfield(X,fieldname,oldvalue{1}); %#ok
    end
end


function	par=wdisp(par,str)

% [par] = wdisp(par,str);
% [par] = wdisp(par);
%			to create a command window <waitbar>
%
% par,str	:	parameter
%--------------------------------------------------------------------------------
%   0,str	:	display <str>	initializes proc
%  >0,str	:	display <str>	accumulates <str>-length+1
%  >0		:	erase	<str>s
%
% note		:
%		<str> must be a CHAR array (string)
%
% example	:
%		t='-\|/';
%	for	i=1:40
%		tl=wdisp( 0,sprintf('time    %c %s',t(rem(i,4)+1),datestr(now)));
%		tl=wdisp(tl,sprintf('count %3d %s',i,repmat('.',1,i)));
%		pause(.1);
%		tl=wdisp(tl);
%	end

% created:
%	us	12-Aug-2000
% modified:
%	us	09-Mar-2002 02:01:37	/ CSSM

	if	~nargin && ~nargout
		help wdisp;
	elseif	nargin == 1
		disp(repmat(char(8),1,par+1));
		par=0;
	elseif	nargin == 2
		disp(str);
		par=par+length(str)+1;
	else
		par=0;
	end
		return;



