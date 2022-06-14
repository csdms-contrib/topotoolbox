function ttscmposter

DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
FD  = FLOWobj(DEM);
S   = STREAMobj(FD,'minarea',100);

FDM = FLOWobj(fillsinks(DEM),'multi');

A   = flowacc(FD);
y = smooth(S,S.y,'k',500,'nstribs',true);
x = smooth(S,S.x,'k',500,'nstribs',true);
k = curvature(S,x,y);
K   = GRIDobj(DEM);
K.Z(:,:) = -inf;
K.Z(S.IXgrid) = k;
K   = dilate(K,ones(3));
K.Z(isinf(K.Z)) = nan;
[~,K] = prcclip(K,10,true);

k = ksn(S,DEM,A);
k = smooth(S,k,'K',200);

DEMf = DEM;
h = fspecial('gaussian',[101 101],31);
DEMf.Z = imfilter(DEMf.Z,h,'symmetric','same');
C1  = curvature(DEMf,'planc');
C2  = curvature(DEMf,'profc');
C3  = curvature(DEMf,'meanc');

c   = chitransform(S,A);
C   = GRIDobj(DEM);
C.Z(S.IXgrid) = c;
C   = dilate(C,ones(5));
C.Z(C.Z==0) = nan;

[U,V] = flowvec(FD);
U.Z = cart2pol(U.Z,V.Z);


cmaps = ttscm;
cmaps = cmaps';
ixx   = [2:4];
cmaps(1,ixx) = {flowdistance(FD) 'flowdistance (up)' []};
cmaps(2,ixx) = {vertdistance2stream(FD,S,DEM) 'vertdistancetostream' []};
cmaps(3,ixx) = {excesstopography(DEM,'maxgradient',20,'kernelsize',31) 'excesstopography' []};
cmaps(4,ixx) = {filter(DEM,'sobel') 'filter (sobel)' [0 200]};
cmaps(5,ixx) = {roughness(filter(DEM)) 'roughness' [0 20]};
cmaps(6,ixx) = {mapfromnal(FD,S,k) 'mapfromnal' []};
cmaps(7,ixx) = {acv(filter(DEM)) 'acv' [2 8]};
cmaps(8,ixx) = {-dilate(sqrt(A),ones(7)) 'flowacc' []};
cmaps(9,ixx) = {dilate(flowdistance(FD,'downstream'),ones(7)) 'flowdistance (down)' []};
cmaps(10,ixx) = {filter(flowconvergence(FDM)) 'flowconvergence' []};
cmaps(11,ixx) = {localtopography(DEM) 'localtopography' []};
cmaps(12,ixx) = {gradient8(DEM) 'gradient8' [0 1]};
cmaps(13,ixx) = {castshadow(DEM,135,10) 'castshadow' [0 1.5]};
cmaps(14,ixx) = {C1 'curvature (plan)' prcclip(C1,1,true)};
cmaps(15,ixx) = {C2 'curvature (profile)' prcclip(C2,1,true)};
cmaps(16,ixx) = {C3 'curvature (mean)' prcclip(C3,1,true)};
cmaps(17,ixx) = {DEM-mean(DEM.Z(:)) 'DEM-mean(DEM)' []};
cmaps(18,ixx) = {distance(DEM,S) 'distance' []};
cmaps(19,ixx) = {K 'stream curvature' []};
cmaps(20,ixx) = {filter(DEM,'wiener') 'filter (wiener)' []};
cmaps(21,ixx) = {shufflelabel(drainagebasins(FD)) 'drainagebasins' []};
cmaps(22,ixx) = {log(flowacc(FDM)) 'log(flowacc) (multi)' []};
cmaps(23,ixx) = {C 'chitransform' []};
cmaps(24,ixx) = {DEM 'DEM' []};
cmaps(25,ixx) = {U 'flowvec' []};
cmaps(26,ixx) = {upslopestats(FD,DEM,'max') 'upslopestats' []};
cmaps(27,ixx) = {aspect(DEM) 'aspect' [0 360]};
cmaps(28,ixx) = {DEM 'anoxia' []};

ixcrop = [199332 598633];
tiledlayout('flow','tilespacing','compact','padding','compact')
for r = 1:size(cmaps,1)
    nexttile
    imageschs(crop(DEM,ixcrop),crop(cmaps{r,2},ixcrop),'colormap', ttscm(cmaps{r,1}),...
        'colorbar',false,'ticklabels','none','caxis',cmaps{r,4},'nancolor',[.4 .4 .4]);
    hold on
    xl = xlim;
    yl = ylim;
    text(xl(1),yl(2),[' ' cmaps{r,1} ' // ' cmaps{r,3} ' '],'FontSize',6,...
        'BackgroundColor','w','Margin',1,'Verticalalignment','top','Clipping','on');
end