function publishtthelp2html

%PUBLISHTTHELP2HTML Publish TopoToolbox mlx help files to html
%
% Syntax
%
%     publishtthelp2html
%
% Description
%
%     HTML files for the documentation of TopoToolbox are published from
%     mlx-files that reside in the mlxfiles-directory. This function
%     publishes and overwrites all HTML-files.
%
% See also: publish
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. November, 2017

% determine the location of this function
fileloc = which('publishtthelp2html.m');
fileloc = fileparts(fileloc);

% set the working directory to those where the mlx-files are in
oldFolder = cd([fileloc filesep 'mlxfiles']);

% list all mlx files
files = dir('*.mlx');

% h = waitbar(0,'Please wait');

% and publish them
nrfiles = numel(files);
parfor r = 1:nrfiles
    
%     waitbar(r/nrfiles,h);
    
    [~,fn] = fileparts(files(r).name);
    % Publish doesn't work with live scripts
    % mydoc{r} = publish(fullfile(files(r).folder,files(r).name),options);
    % Thanks to 
    % https://de.mathworks.com/matlabcentral/answers/282820-programmatically-run-and-export-live-script
    % following undocumented approach works
    
    % matlab.internal.liveeditor.executeAndSave(files(r).name);
    matlab.internal.liveeditor.openAndConvert(files(r).name, [fileloc filesep fn '.html']);
end
% close(h);
% Reset working directory back to previous working directory
builddocsearchdb(fileloc)
cd(oldFolder);