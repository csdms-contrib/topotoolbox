function C = ttbenchmark


load exampledem
DEM = GRIDobj(X,Y,dem);
clear X Y dem



for r = 1:6;
    
    C{r+1,1} = DEM.size;
    C{r+1,2} = prod(C{r+1,1});
    
    C{1,3} = 'FLOWobj';
    F = @() FLOWobj(DEM,'PreProcess','Carve',...
                       'mex',false);    
    C{r+1,3} = timeit(F);
    
    C{1,4} = 'FLOWobj (mex)';
    F = @() FLOWobj(DEM,'PreProcess','Carve',...
                       'mex',true);    
    C{r+1,4} = timeit(F);
    
    C{1,5} = 'flowacc';
    FD = FLOWobj(DEM,'PreProcess','Carve',...
                       'mex',true); 
    F = @() flowacc(FD);
    C{r+1,5} = timeit(F);
    
    
    % resample dem
    
    DEM = resample(DEM,2);    
end


loglog(cell2mat(C(2:end,2)),cell2mat(C(2:end,3)),'-sk')
hold on
loglog(cell2mat(C(2:end,2)),cell2mat(C(2:end,4)),'-sr')

xlabel('# of cells');
ylabel('time [s]');



