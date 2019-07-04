%TODO define seperate BCs for edges (eg. combination of Neuman and Dirichlet etc. )

function u_g= BC_Disturbance(u_g,distSites,BC)

if any(strcmp(distSites,{'l','lr','tl','bl','tblr'}))
    u_g(BC.BC_indices.leftRow)=u_g(BC.BC_indices.leftRow).*rand(size(u_g(BC.BC_indices.leftRow)))*BC.BC_dir_Dist_Value;
end
if any(strcmp(distSites,{'r','lr','tr','br','tblr'}))
    u_g(BC.BC_indices.rightRow)=u_g(BC.BC_indices.rightRow).*rand(size(u_g(BC.BC_indices.rightRow)))*BC.BC_dir_Dist_Value;
end
if any(strcmp(distSites,{'t','tr','tl','bt','tblr'}))
    u_g(BC.BC_indices.topRow)=u_g(BC.BC_indices.topRow).*rand(size(u_g(BC.BC_indices.topRow)))*BC.BC_dir_Dist_Value;
end
if any(strcmp(distSites,{'b','br','tb','bl','tblr'}))
    u_g(BC.BC_indices.botRow)=u_g(BC.BC_indices.botRow).*rand(size(u_g(BC.BC_indices.botRow)))*BC.BC_dir_Dist_Value; 
end

end