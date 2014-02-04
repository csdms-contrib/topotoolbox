function S = extractconncomps(S,cc)

% interactive stream network selection
%
% Syntax
%
%     S2 = extractconncomps(S)
%
% Description
%
%     extractconncomps displays a figure and plots the stream network S.
%     Numbers at the stream outlets refer to single trees in the stream
%     network and the user can choose which should be retained in the
%     output STREAMobj S2 by entering the numbers in the dialog box.
%
% Input arguments
%
%     S    instance of a STREAMobj
%
% Output arguments
%
%     S2   instance of a STREAMobj
%
% See also: STREAMobj/klargestconncomps
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013 

nrc = numel(S.x);
M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);

[~,p,~,r] = dmperm(M | M' | speye(nrc));
nc = length(r) - 1;

% label matrix
L = zeros(nrc,1);
for tt = 1:nc;
    L(p(r(tt):r(tt+1)-1)) = tt;
end

if nargin == 1;

hFig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
hAx  = axes('Parent',hFig);
cmap = jet(nc);

for r = 1:nc;
    M2 = spdiags(L==r,0,nrc,nrc)*double(M);
    [x,y] = gplot(M2,[S.x S.y]);
    outlet = find((sum(M2,1)' > 0) & (sum(M2,2) == 0));
    
    plot(hAx,x,y,'Color',cmap(r,:));
    hold on
    plot(hAx,S.x(outlet),S.y(outlet),'*','Color',cmap(r,:));
    text(S.x(outlet),S.y(outlet),[' ' num2str(r)]);
end

axis image
set(hAx,'Xtick',[],'Ytick',[]);

prompt = {'Enter comma separated list of component indices to extract:'};

dlg_title = 'Extract connected components';
num_lines = 1;
options.Resize='on';
answer = inputdlg(prompt,dlg_title,num_lines,{''},options);

% evaluate answer
c = regexp(answer,'(?:\d*)','match');
c = str2double(c{1});
else
    c = cc;
    if c > nc;
        error('TopoToolbox:wronginput','The supplied number exceeds the number of components');
    end
end

% adapt new STREAMobj to the reduced network
L     = ismember(L,c);
I     = L(S.ix);
S.ix  = S.ix(I);
S.ixc = S.ixc(I);

IX    = cumsum(L);
S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

S.x   = S.x(L);
S.y   = S.y(L);
S.IXgrid   = S.IXgrid(L);






    



