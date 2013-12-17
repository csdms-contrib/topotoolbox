function C = ttbenchmark

% 

load exampledem
DEM = GRIDobj(X,Y,dem);
clear X Y dem



for r = 1:6;
    
    C{1,1} = 'rows';
    C{r+1,1} = DEM.size(1);
    C{1,2} = 'columns';
    C{r+1,2} = DEM.size(2);
    C{1,3} = 'element count';
    C{r+1,3} = prod(DEM.size);
    
    C{1,4} = 'FLOWobj [s]';
    F = @() FLOWobj(DEM,'PreProcess','Carve',...
                       'mex',false);    
    C{r+1,4} = timeit(F);
    
    C{1,5} = 'FLOWobj (mex) [s]';
    F = @() FLOWobj(DEM,'PreProcess','Carve',...
                       'mex',true);    
    C{r+1,5} = timeit(F);
    
    C{1,6} = 'flowacc (FLOWobj) [s]';
    FD = FLOWobj(DEM,'PreProcess','Carve',...
                       'mex',true); 
    F = @() flowacc(FD);
    C{r+1,6} = timeit(F);
    
    C{1,7} = 'flowacc (M) [s]';
    M = FLOWobj2M(FD);
    F = @() flowacc(M,DEM.size);
    C{r+1,7} = timeit(F);
    
    C{1,8} = 'speed increase [%]';
    C{r+1,8} = C{r+1,7}/C{r+1,6};
    
    % resample dem   
    DEM = resample(DEM,2);    
end


loglog(cell2mat(C(2:end,2)),cell2mat(C(2:end,3)),'-sk')
hold on
loglog(cell2mat(C(2:end,2)),cell2mat(C(2:end,4)),'-sr')

xlabel('element count');
ylabel('time [s]');



