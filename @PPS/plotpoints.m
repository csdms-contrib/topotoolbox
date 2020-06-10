function h = plotpoints(P,varargin)

%PLOTPOINTS plot points of PPS
%
% Syntax
%
%     plotpoints(P)
%     plotpoints(P,'pn',pv,...)
%     h = plotpoints(...)
%
% Description
%
%     plotpoints plots the points of an instance of PPS. The function uses
%     the build-in function scatter. Thus, it accepts all parameter
%     name/value pairs that are accepted by scatter.
%     
% Input arguments
%
%     P    instance of PPS
%
%     Parameter name/value pairs (see scatter)
%
%     'Marker'     Marker type (see scatter)
%     'MarkerEdgeColor'  Marker edge color (see scatter) (default = 'k')
%     'MarkerFaceColor'  Marker face color (see scatter) (default = 'w')
%     'LineWidth'  Width of the marker edge
%     'SizeData'   Default=20. If scalar, then markers will be plotted
%                  using this size (see scatter). If vector, then marker 
%                  size will vary according to the values in this vector. 
%                  Marker size will be automatically scale to vary within
%                  the bounds set by 'MinSize' and 'MaxSize'.
%     'MinSize'    see 'SizeData' (default = 5)
%     'MaxSize'    see 'SizeData' (default = 75)
%     'MarkerFaceAlpha'  Transparency (default = 0.5). (see scatter)
%     'covariate'  Valid covariate (e.g. GRIDobj, node-attribute list) that
%                  will be used for coloring the points
%     'marks'      Vector with marks that will be used for coloring the
%                  points. 'covariate' and 'marks' cannot be set jointly.
%
% Example
%
%
%
%
% See also: PPS/plot, PPS/plotdz, STREAMobj/plot, STREAMobj/plotc
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'Marker','o');
addParameter(p,'MarkerEdgeColor','k');
addParameter(p,'MarkerFaceColor','w');
addParameter(p,'MarkerEdgeAlpha',0.5);
addParameter(p,'MarkerFaceAlpha',0.5);

addParameter(p,'LineWidth',1)
addParameter(p,'SizeData',20);
addParameter(p,'MinSize',5);
addParameter(p,'MaxSize',75);
addParameter(p,'covariate',[]);
addParameter(p,'marks',[]);

% Parse
parse(p,varargin{:});

Results = p.Results;
if ~isempty(Results.covariate) && ~isempty(Results.marks)
    error('Cannot set coloring to both covariate and marks')
end

if ~isempty(Results.covariate)
    c = getmarks(P,Results.covariate);
    Results = rmfield(Results,'MarkerFaceColor');
    Results.CData = +c;
    
elseif ~isempty(Results.marks)
    c = Results.marks;
    Results = rmfield(Results,'MarkerFaceColor');
    Results.CData = +c;
end

if isscalar(Results.SizeData)
    % all good. All points have the same size
else
    if numel(Results.SizeData) == npoints(P)
    else
        Results.SizeData = getmarks(P,Results.SizeData);
    end
    % scale between 1 and 20
    Results.SizeData = Results.SizeData - min(Results.SizeData);
    Results.SizeData = Results.SizeData./max(Results.SizeData);
    Results.SizeData = Results.SizeData*(Results.MaxSize - Results.MinSize) + Results.MinSize;
end
        
    
Results = rmfield(Results,{'covariate','marks','MaxSize','MinSize'});

xy = P.ppxy;

pnpv = expandstruct(Results);
% pn = fieldnames(Results);
% pv = struct2cell(Results);
% 
% pnpv = [pn pv];
% pnpv = pnpv';
% pnpv = pnpv(:)';

ht = scatter(xy(:,1),xy(:,2),'filled',pnpv{:});

if nargout == 1
    h = ht;
end

end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end