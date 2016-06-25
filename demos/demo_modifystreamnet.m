function demo_modifystreamnet

% demo on geometric modifications of stream networks
%
% Syntax
%
%     demo_modifystreamnet
%
% Description
%
%     This demo shows most of the available, automated tools to modify the 
%     planform geometry of stream networks in TopoToolbox. The function
%     requires no input and output arguments.
%
% See also: STREAMobj/modify, STREAMobj/trunk, STREAMobj/klargestconncomps
%           STREAMobj/removeshortstreams
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. April, 2015

narginchk(0,0)
nargoutchk(0,0)

load exampledem
DEM = GRIDobj(X,Y,dem);
FD  = FLOWobj(DEM,'preprocess','carve');
S   = STREAMobj(FD,'minarea',100);
c1   = [.5 .5 .5];
c2   = 'r';
lw   = 1.5;

%% klargestconncomps
subplot(3,5,1);
plot(S,'color',c1);
hold on
S2 = klargestconncomps(S,1);
plot(S2,c2,'LineWidth',lw);
setaxisright
title('klargestconncomps(S,1)');
hold off

%% removeshortstreams(S,1000)
subplot(3,5,2);
plot(S,'color',c1);
hold on
S2 = removeshortstreams(S,1000);
plot(S2,c2,'LineWidth',lw);
setaxisright
title('removeshortstreams(S,1000)');
hold off

%% trunk
subplot(3,5,3);
plot(S,'color',c1);
hold on
S2 = trunk(S);
plot(S2,c2,'LineWidth',lw);
setaxisright
title('trunk(S)');
hold off

%% modify(S,'distance',10000)
subplot(3,5,4);
plot(S,'color',c1);
hold on
S2 = modify(S,'distance',10000);
plot(S2,c2,'LineWidth',lw);
setaxisright
title('modify(S,''distance'',10000)');
hold off

%% modify(S,'distance',[5000 10000])
subplot(3,5,5);
plot(S,'color',c1);
hold on
S2 = modify(S,'distance',[5000 10000]);
plot(S2,c2,'LineWidth',lw);
setaxisright
title({'modify(S,''distance'',...';'[5000 10000])'});
hold off

%% modify(S,'streamorder',1)
subplot(3,5,6);
plot(S,'color',c1);
hold on
S2 = modify(S,'streamorder',1);
plot(S2,c2,'LineWidth',lw);
setaxisright
title('modify(S,''streamorder'',1)');
hold off

%% modify(S,'streamorder','<=2')
subplot(3,5,7);
plot(S,'color',c1);
hold on
S2 = modify(S,'streamorder','<=2');
plot(S2,c2,'LineWidth',lw);
setaxisright
title({'modify(S,''streamorder'',...';'''<=2'')'});
hold off

%% modify(S,'upstreamto',DEM>=1600)
subplot(3,5,8);
plot(S,'color',c1);
hold on
S2 = modify(S,'upstreamto',DEM>=1600);
plot(S2,c2,'LineWidth',lw);
setaxisright
title({'modify(S,''upstreamto'',...';'DEM>=1600)'});
hold off

%% modify(S,'downstreamto',DEM<=1600)
subplot(3,5,9);
plot(S,'color',c1);
hold on
S2 = modify(S,'downstreamto',DEM<=1600);
plot(S2,c2,'LineWidth',lw);
setaxisright
title({'modify(S,''downstreamto'',...';'DEM<=1600)'});
hold off

%% modify(S,'tributaryto',trunk(S))
subplot(3,5,10);
plot(S,'color',c1);
hold on
S2 = modify(S,'tributaryto',trunk(S));
plot(S2,c2,'LineWidth',lw);
setaxisright
title({'modify(S,''tributaryto'',...';'trunk(S))'});
hold off

end

function setaxisright(ax)
if nargin == 0;
    ax = gca;
end

axis(ax,'image');
set(ax,'Xticklabel',[]);
set(ax,'Yticklabel',[]);
end
