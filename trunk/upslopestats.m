function S = upslopestats(M,VAR,type)

% upslope statistics of a variable based on the flow direction matrix
%
% Syntax
%
%     S = upslopestats(M,VAR,type)
%
% Description
%
%     upslopestats returns various statistical parameters such as mean and
%     standard deviation of each cell's upslope located samples VAR based on
%     the flow direction matrix M. Note that while mean and standard 
%     deviation are calculated fairly fast, the calculation of minimum and 
%     maximum is non-vectorized and usually takes longer to evaluate.
%
% Input
%
%     M         sparse flow direction matrix (single and multi)
%     VAR       variable (e.g. slope calculated by gradient8)
%     type      'mean' (default), 'std' standard deviation, 
%               'var' variance, 'min', 'max'               
%
% Output
% 
%     S         array same size as VAR
%
% Example
%
%     % calculate the network wetness index (Lane et al. 2004) 
%     load exampleDEM
%     M = flowdir_single(dem);
%     A = flowacc(M,size(dem));
%     G = gradient8(dem,90);
%     % topographic wetness index
%     T = twi(A,G);
%     % identify channels as a simple threshold of 
%     % drainage area
%     CHANNELS = A>500;
%     % set the TWI of channels to inf
%     T(CHANNELS) = inf;
%     % upstream minima imposition (note that M in this
%     % case is transposed!)
%     Tn = upslopestats(M',T,'min');
%
%
% References
%
%     Lane, S. N., Brookes, C. J., Kirkby, M. J. & Holden, J. A, 2004: A
%     network-index-based version of TOPMODEL for use with high-resolution 
%     digital topographic data Hydrological Processes, 18, 191-201. 
%
% See also: FLOWDIR, FLOWACC
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 23. May, 2011


if nargin == 2;
    type = 'mean';
else
    type = lower(type);
end





switch lower(type)
    case 'sum'
        S = flowacc(M,VAR);
    case 'mean'
        S = flowacc(M,VAR)./flowacc(M,size(VAR));
    case {'std','var'}
        n = flowacc(M,size(VAR));
        % sum
        S = flowacc(M,VAR);
        % variance
        S = 1./(n-1) .* ((flowacc(M,VAR.^2) - 1./n.*(S.^2)));
        S(n==1) = 0;
        
        switch lower(type)
            case 'std'
                S = sqrt(S);
        end

    case {'min','max'}
        
        % get indices        
        [ic,icd]   = find(M);
        % get non-receivers
        nonrec = full((sum(M,1) == 0)');
        
        switch lower(type)
            case 'min'
                relop = @lt;
                
                t = sortrows([nonrec(ic) ic icd VAR(ic) VAR(icd)],[-1 4 5]);
            case 'max';
                relop = @gt;
                
                t = sortrows([nonrec(ic) ic icd VAR(ic) VAR(icd)],[-1 -4 -5]);
        end
        
        nonrec = ~t(:,1);
        ic     = t(:,2);
        icd    = t(:,3);
        v      = t(:,4);
        
        [~,icdd] = ismember(icd,ic);
        
        
        % initialize stat variable
        S  = VAR;
        S(ic(nonrec)) = v(nonrec);
        
        NG  = find(nonrec);
        % nr of nonreceivers
        nNG = numel(NG);
        
        NGix = 1;     
        IX   = 1;
        
        
        % not vectorized
        while NGix <= nNG;
            
            % distance
            ln = S(ic(IX));
            
            if relop(ln,S(icd(IX)));
                S(icd(IX)) = ln;
                
                % next row to work on
                IX = icdd(IX);
                
                flaggoon = IX~= 0;
            else
                flaggoon = false;
            end
            
            if ~flaggoon;
                
                NGix = NGix+1;
                
                if NGix>nNG;
                else
                    IX  = NG(NGix);
                end
            end
        end
        
        
end

    

    
