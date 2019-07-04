function cmap = loadcmap(cmaptype)

%LOADCMAP load colormap

p = fileparts(mfilename('fullpath'));
cmap = load([p filesep cmaptype filesep cmaptype '.mat']);
cmap = cmap.(cmaptype);