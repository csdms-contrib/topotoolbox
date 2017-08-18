function plotsegmentgeometry(S,segment)

%PLOTSEGMENTGEOMETRY Plot segment geometry obtained from the function networksegment
%
% Syntax
%
%     plotsegmentgeometry(S,segment)
%
% Description
%     
%     This function plot the river segments orientation, euclidian length
%     flow length, sinuosity, segment elevation difference, slope, drainage
%     area, mean steepness (assuming m/n=0.5) and strahler order
%
% Input arguments
%
%     S         stream network (STREAMobj)  
%     segment   River segments object
%
% Example (see function networksegment)
%
% 
% See also: STREAMobj, networksegment
% 
% Author: Philippe Steer (philippe.steer[at]univ-rennes1.fr)
% Date: 11. May, 2016

%% Plot
% Number of colors
nhs = 256;
lwidth=2;

% Plot Strahler order
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(1,nanmax(segment.strahler),nhs);plot(S.x,S.y,'k.','MarkerSize',1);hold on;title('Strahler order');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.strahler(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',2);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot Orientation
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,180,nhs);
subplot(1,4,[1 3]);plot(S,'k');hold on;title('Orientation');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.angle(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;caxis([cmap_vec(1) cmap_vec(end)]);
subplot(1,4,4);
thetas=(0:1:180)'*pi/180; r = 1; [x,y]=pol2cart(thetas,r);  
for i=1:numel(x)-1;[~,indcmap]=min(abs(cmap_vec-thetas(i).*180/pi));X=[0 x(i) x(i+1) 0];Y=[0 y(i) y(i+1) 0];patch(X,Y,cmap(indcmap,:),'EdgeAlpha',0);end; axis equal off;

% Plot Euclidian length
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,nanmax(segment.length),nhs);plot(S,'k');hold on;title('Euclidian length');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.length(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot flow length
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,nanmax(segment.flength),nhs);plot(S,'k');hold on;title('Flow length');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.flength(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot sinuosity
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(1,nanmax(segment.sinuosity)-1,nhs);plot(S,'k');hold on;title('Sinuosity');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.sinuosity(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot segment slope
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,nanmax(segment.slope),nhs);plot(S,'k');hold on;title('Segment slope');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.slope(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot slope along the flow path
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,nanmax(segment.fslope),nhs);plot(S,'k');hold on;title('Flow path slope');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-segment.fslope(i)));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot mean drainage area
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,log10(nanmax(segment.Amean)),nhs);plot(S,'k');hold on;title('Mean drainage area - log_{10}');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-log10(segment.Amean(i))));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

% Plot mean steepness
figure;cmap=colormap(jet(nhs));cmap_vec=linspace(0,log10(nanmax(segment.ksn)),nhs);plot(S,'k');hold on;title('Mean steepness - log_{10}');
for i=1:segment.n; [~,indcmap]=min(abs(cmap_vec-log10(segment.ksn(i))));plot([S.x(segment.ix(i,1)) S.x(segment.ix(i,2))],[S.y(segment.ix(i,1)) S.y(segment.ix(i,2))],'-','Color',cmap(indcmap,:),'LineWidth',lwidth);end;colormap(cmap);axis equal tight;colorbar;caxis([cmap_vec(1) cmap_vec(end)]);

end