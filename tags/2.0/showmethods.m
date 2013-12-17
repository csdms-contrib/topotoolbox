function showmethods(classname)

% displays class method names and H1 lines in the command line
%
% Syntax
%
%     showmethods(classname)
%
% Description
%
%     showmethods outputs a list of methods available for a specific class
%     (e.g. GRIDobj) in the command window. Only the methods are listed
%     that feature an H1 line, e.g., the first line of the help text block
%     in a function.
%
% Input arguments
%
%     classname     string with class name (e.g. GRIDobj, FLOWobj)
%
% Example
%
%     showmethods('GRIDobj')
%
% See also: methods, properties, class
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. June, 2013

narginchk(1,1)

m = methods(classname);

C = cellfun(@(x) numel(x),m,'uniformoutput',true);

maxmethcharacter = max(C);

for r = 1:numel(m)
    s = which([classname '/' m{r}]);
    fileID = fopen(s,'r');
    H1line = false;
    
    try
        while ~H1line && ~feof(fileID)
            tline = fgetl(fileID);
            if isempty(tline)
                continue
            end
            if strcmp(tline(1),'%');
                H1line = true;
            end
        end
        
        if feof(fileID)
            fclose(fileID);
            continue
        end
        
        methodstr = m{r};
        if numel(methodstr)<= maxmethcharacter;
            methodstr = [methodstr repmat(' ',1,maxmethcharacter - numel(methodstr))]; %#ok<AGROW>
        end
        
        h1str = tline(3:end);
                        
        disp([ upper(methodstr) ' : ' h1str]);
        fclose(fileID);
        
    catch
        disp(m{r})
    end
    
    
end
            
           


