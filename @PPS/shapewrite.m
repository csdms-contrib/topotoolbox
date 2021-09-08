function shapewrite(P,filename,varargin)

%SHAPEWRITE write point pattern to shapefile
%
% Syntax
%
%     shapewrite(P,filename)
%     shapewrite(P,filename,pn,pv,...)
%
% Description
%
%     shapewrite writes both the point pattern and the stream network to a
%     shapefile. The points are written to a point shapefile
%     and the stream network to a line shapefile. 
%
% Input arguments
%
%     P         point pattern (PPS)
%     filename  filename of the output shapefile (no *.shp required)
%
%     Parameter name/value pairs
%
%     'streams'  {false} or true. If false, shapewrite will not write the 
%                shapefile of the stream network
%     'postfix'  {'_stream'} postfix for the stream network
%     'marks'    attribute data for each point. This can be either a cell
%                array of size 2xn where n is the number of attributes. The
%                elements in the cell array must contain the name of the
%                attribute and the data. Data can be GRIDobj or node
%                attribute lists. Alternatively, one can supply a table
%                with marks.
%
%
% See also: STREAMobj/STREAMobj2mapstruct 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

[filepath,name,~] = fileparts(filename);

p = inputParser;
p.FunctionName = 'PPS/shapewrite';
addParameter(p,'streams',false);
addParameter(p,'postfix','_stream');
addParameter(p,'marks',[]);
parse(p,varargin{:});


xy = P.ppxy;
id = (1:npoints(P))';
sh = struct('Geometry','Point',...
            'X',num2cell(xy(:,1)),...
            'Y',num2cell(xy(:,2)),...
            'id',num2cell(id));
if ~isempty(p.Results.marks)
    if iscell(p.Results.marks)
        for r = 1:2:numel(p.Results.marks)
            mark = num2cell(p.Results.marks{r+1});
            [sh.(p.Results.marks{r})] = mark{:};
        end
    elseif istable(p.Results.marks)
        sh = struct2table(sh);
        sh = [sh p.Results.marks];
        sh = table2struct(sh);
%         t = struct2table(p.Results.marks);
%         sh = [sh;t];
    end
end

shapewrite(sh,fullfile(filepath,name))

if p.Results.streams
    MS = STREAMobj2mapstruct(P.S);
    shapewrite(MS,fullfile(filepath,[name p.Results.postfix]));
end
