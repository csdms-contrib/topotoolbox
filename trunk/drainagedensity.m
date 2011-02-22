function [d,a,l,D] = drainagedensity(X,Y,Mstreams,DB)

% calculate drainage density of individual drainage basins
%
% Syntax
%
%     [d,a,l] = drainagedensity(X,Y,Mstreams,DB)
%
% Description
%
%     Drainage density d is defined by d=sum(L)/Ad where L is the stream 
%     length and Ad is the drainage basin area. Drainage density is 
%     regarded as a fundamental measure of drainage network geometry. It 
%     describes the basin dissection by surface streams (Knighton 1998).
%
% Input
%
%     X,Y       coordinate matrices (as produced by meshgrid)
%     Mstreams  single flow direction matrix for streams only (third output
%               argument of streamorder)
%     DB        drainage basins (as returned by drainagebasins) or any
%               label matrix.
% 
% Output
%
%     d         drainage density for each drainage basin
%     a         drainage basin area
%     l         sum of stream length in each basin
%     DENS      map of drainage densities
%
% Example 
%
%     load exampleDEM
%     % calculate flow accumulation and direction
%     [A,M] = ezflowacc(X,Y,dem,'type','single');
%     % let's simply assume that channels start where
%     % A is larger than 100;
%     W = A>100;
%     % and calculate the strahler stream order.
%     [S,nodes,Mstreams] = streamorder(M,W);
%     % then calcuate the extent of drainage basins of 2. order
%     DB = drainagebasins(M,S,2);
%     % and then calcuate drainage density for each basin
%     [d,a,l,DENS] = drainagedensity(X,Y,Mstreams,DB);
%     hist(d)
%     xlabel('drainage density [m/m^2]')
%     ylabel('absolute frequency')
%     figure
%     imagesc(X(1,:),Y(:,2),DENS); axis image; axis xy
%     hold on
%     gplot(Mstreams,[X(:) Y(:)],'k');
%     hold off
%     colorbar
%     title('drainage density')
%     
%
%
%
% See also: FLOWDIR_SINGLE, STREAMORDER, REGIONPROPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 10. July, 2009


if ~isequal(size(DB),size(X),size(Y))
    error('TopoToolbox:incorrectinput',...
          'DB and W must have same size')
end

DB(isnan(DB)) = 0;
cs = abs(Y(1)-Y(2));

% if someone might have a different coordinate system
if cs == 0
    cs = abs(X(1)-X(2));
end

% calculate stream distances
[ic,icd] = find(Mstreams);
D        = zeros(size(X));
D(ic)    = hypot(X(ic)-X(icd),Y(ic)-Y(icd));


% check image processing toolbox version
if verLessThan('images','6.3')
    s = regionprops(DB,'area');
    a = [s.Area];
    a = a(:);
    DB = DB(:);
    I  = DB>0;
    l = accumarray(DB(I),D(I));
       
else
    s = regionprops(DB,D,'area','Pixelvalues');
    a = [s.Area];
    a = a(:);
    l = cell2mat(cellfun(@(x) sum(x),{s.PixelValues}','Uniformoutput',false));
    
end

% calculate density
d = l./(a*(cs^2));

% produce forth output if wanted
if nargout==4
    if ~exist('I','var')
    	  I = DB>0;
    end
    D    = zeros(size(X));
    D(I) = d(DB(I));
    D(~I)= nan;
end





