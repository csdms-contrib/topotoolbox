function br = sidebranching(S,outputformat)

%SIDEBRANCHING side branching classification according to Tokunaga (1978)
%
% Syntax
%
%     br = sidebranching(S)
%     br = sidebranching(S,outputformat)
%
% Description
%
%     The Strahler stream order is a common way to classify rivers.
%     However, occasionally this information may not suffice to describe
%     the branching behavior of river networks. Instead, we may wish to
%     quantify which river-orders intersect how often with other river
%     orders. E.g., if a first-order stream intersects with another
%     first-order stream, we may wish to denote that stream as '11'; if a
%     first-order stream intersects with a third-order stream, we denote
%     that one as '13'. This is the method of Tokunaga (1978) which extends
%     on the Strahler (1957) order system (Pelletier and Turcotte, 2000).
%
% Input arguments
%
%     S       STREAMobj
%     outputformat   output format can be 'nalstr' (node attribute list that
%             comes as cell array with strings), 'nalnum' (node attribute
%             list with classifications as double), or 'mapstruct'.
%
% Output arguments
%
%     br      see outputformat
%
% Example
%
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     MS = sidebranching(S,'mapstruct');
% 
%     % Visualize results
%     [unbr,~,loc] = unique({MS.sidebranch});
%     c = lines(numel(unbr));
%     hold on
%     for r = 1:numel(MS);
%         h(loc(r)) = plot(MS(r).X,MS(r).Y,'Color',c(loc(r),:),'LineWidth',MS(r).streamorder);
%     end
%     legend(h,unbr,'location','best')
%
% See also description here: 
%
%     https://topotoolbox.wordpress.com/2015/07/21/side-branching/
% 
% References
%
%     Pelletier, J.D., Turcotte, D.L. 2000: Shapes of river networks and
%     leaves: are they statistically similar? Phil. Trans. R. Soc. Lond. B,
%     355, 307-311.
%
%     Strahler, A.N. 1957: Quantitative analysis of watershed
%     geomorphology. Am. Geophys. Union Trans. 38, 913-920.
%
%     Tokunaga, E. 1978: Consideration on the composition of drainage
%     networks and their evolution. Geogr. Repl. Tokyo Metro. Univ. 13,
%     1-27.
%
% See also: STREAMobj2mapstruct, STREAMobj/streamorder
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. July, 2015


% check input arguments
if nargin==1;
    outputformat = 'nalstr';
else
    outputformat = validatestring(outputformat,{'nalstr' 'nalnum' 'mapstruct'});
end

% get confluences
confl = streampoi(S,'confluences','logical') | streampoi(S,'outlets','logical');

% calculate stream order but record at the same the branching at
% confluences
s  = +streampoi(S,'Channelheads','logical');
br = cell(size(s));
br = cellfun(@(x) '', br,'UniformOutput',false);

% calculate stream order and record type of confluences one stream pixel
% upstream of confluences
for r = 1:numel(S.ix);
    if s(S.ixc(r)) == 0 || s(S.ixc(r)) <  s(S.ix(r));
        s(S.ixc(r)) = s(S.ix(r));
    elseif s(S.ixc(r)) ==  s(S.ix(r));
        s(S.ixc(r)) = s(S.ix(r))+1;
    end
 
    if confl(S.ixc(r));
        br{S.ixc(r)} = [num2str(s(S.ix(r))) br{S.ixc(r)}];
    end
end

% sort side branching ascending
br(confl) = cellfun(@(x) sort(x,'ascend'),br(confl),'UniformOutput',false);

% check side branches and remove those where higher order river join lower
% order rivers
for r = 1:numel(S.ix);
    if confl(S.ixc(r));
        s1 = br{S.ixc(r)};
        s1 = str2double(s1(1));
        if s(S.ix(r)) == s1
            br{S.ix(r)} = br{S.ixc(r)};
        end
    elseif confl(S.ix(r));
        br{S.ix(r)} = '';
    end
end

% follow and record side branching upstream
for r = numel(S.ix):-1:1;
    if isempty(br{S.ix(r)})
        br(S.ix(r)) = br(S.ixc(r));
    end
end

% generate output
switch outputformat
    case 'nalnum'
        br = cellfun(@(x) str2double(x),br,'UniformOutput',true);
    case 'mapstruct'
        MS  = STREAMobj2mapstruct(S);
        xy  = zeros(numel(MS),2);
        for r = 1:numel(MS);
            xy(r,1) = MS(r).X(1);
            xy(r,2) = MS(r).Y(1);
        end
        
        [~,locb] = ismember(xy,[S.x S.y],'rows');
        for r = 1:numel(MS);
            MS(r).sidebranch = br{locb(r)};
        end
        
        br = MS;
end
    
        