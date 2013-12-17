function SH = castshadow(DEM,azid,altd)

% cast shadow
%
% Syntax
%
%     SH = castshadow(DEM,azid,altd)
%
% Description
%     
%     castshadow calculates the shadow created on a form next to a surface 
%     that is turned away from the source of light. 
%
% Input arguments
%
%     DEM    digital elevation model (class: GRIDobj)
%     azid   azimuth angle in degrees of light source (clockwise from top)
%     altd   elevation angle in degrees of light source above horizon
%
% Output arguments
%
%     SH     values between 0 and 1 with values of 1 indicating shadowed 
%            and 0 indicating non-shadowed areas. Values between 0 and 1
%            can be interpreted to indicate a degree of shadowing. 
%
% Example
%     
%     % Sun azimuth and elevation for Zürich, 30. April, 2013
%     % Data taken from http://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
%     
%     % hour
%     h = 6:19;
%     azid = [75.72 86.44 97.65 110.24 125.52 ...
%             145.09 169.68 196.37 219.91 238.38 ...
%             252.92 265.11 276.15 286.87];
% 
%     altd = [7.14 17.1 27.23 37.07 46.05 53.22 ...
%             57.17 56.64 51.83 44.15 34.92 24.98 ...
%             14.87 5.05];
% 
%     % DEM = GRIDobj('D:\WDmatlab\mymfiles\demanalysis\examples\grindelwald.txt');
%     for r = 1:numel(h); 
%         I = castshadow(2.*DEM,azid(r),altd(r)); 
%         imageschs(DEM,I,'colormap',flipud(gray)*.8 + .2,'azimuth',azid(r),'altitude',altd(r));
%         title(num2str(h(r)));
%         pause(.5);
%     end
%
% 
% See also: GRIDobj/imageschs, GRIDobj/hillshade  
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. April, 2013



azid = mod(-azid,360);
% force values in the DEM to be larger than zero
% because imrotate fills loose image ends with zeros
DEM.Z =  - min(DEM.Z(:)) + DEM.Z + 1;

%
% affine transformation matrix
A = [cosd(azid) sind(azid); ...
    -sind(azid) cosd(azid)];
A(3,3) = 1;    
tform = maketform('affine',A);
% image transformation
[demr,xdata,ydata]  = imtransform(DEM.Z,tform,'bilinear',...
    'Udata',[1 DEM.size(2)],'Vdata',[1 DEM.size(1)]);

% calculate height of sunbeams
lowering = DEM.cellsize*tand(altd);
demrcopy = demr;
for r = 2:size(demr,1);
    demr(r,:) = max(demr(r-1,:) - lowering,demr(r,:)).*(demr(r,:)>0);
end
demr = single(demr>demrcopy);

% backtransform image 
azid = -azid;
A = [cosd(azid) sind(azid); ...
    -sind(azid) cosd(azid)];
A(3,3) = 1;
demr = imtransform(demr,maketform('affine',A),'bilinear',...
    'Udata',xdata,'Vdata',ydata,...
    'xdata',[1 DEM.size(2)],'ydata',[1 DEM.size(1)]);

% write to GRIDobj
SH = DEM;
SH.Z = demr;




    