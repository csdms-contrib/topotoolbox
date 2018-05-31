function Z= setBC(Z,BC)
% Function to set Boundary conditions
%
% Syntax
%
%       Z= setBC(Z,BC)
%
% Description
%
%       Function set the boundary conditions of a simulated domain (Z)
%
% Input
%
%       Z         Simualted domain
%       BC        Boundary conditions as an instance of BC
% 
% 
% Output
%
%       Z         Simualted domain with updated boundaries
%
% Example
%
%
% See also:
%
% Authors: Benjamin Campforts (benjamin.campforts[at]ees.kuleuven.be)
%          Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date: 10. June, 2016

switch BC.type
    case 'Periodic'
        switch BC.nbGhost
            case 1
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow    BC.BC_indices.botRow  BC.BC_indices.topRow  ])=Z([ BC.BC_indices.leftRow_2  BC.BC_indices.rightRow_2  BC.BC_indices.topRow_2   BC.BC_indices.botRow_2 ]);
            case 2
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2 ])=Z([ BC.BC_indices.leftRow_2  BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_3  BC.BC_indices.rightRow_3  BC.BC_indices.topRow_2   BC.BC_indices.botRow_2  BC.BC_indices.topRow_3  BC.BC_indices.botRow_3 ]);
        end
    case 'Dirichlet'
        if ~isequal(size(BC.dir_value),[1 1])
            error('p.BC_dir_value should be a single value, otherwise, use the Dirichlet_Matrix BC');
        end
        switch BC.nbGhost
            case 1
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.botRow  BC.BC_indices.topRow ])=BC.dir_value;
            case 2
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2])=BC.dir_value;
        end
    case 'Dirichlet_Matrix'
        if ~isequal(size(BC.dir_value),size(Z))
            error('p.BC_dir_value should be a (sparse) matrix with the same dimensions as the modelled domain (DEM.Z)');
        end
        switch BC.nbGhost
            case 1
                Z([ BC.BC_indices.topRow  BC.BC_indices.botRow  BC.BC_indices.leftRow  BC.BC_indices.rightRow])=BC.dir_value([ BC.BC_indices.topRow  BC.BC_indices.botRow  BC.BC_indices.leftRow  BC.BC_indices.rightRow]);
            case 2
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2])=BC.dir_value([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2]);
        end
    case 'Dirichlet_Matrix_Ini'
        if ~isequal(size(BC.dir_value),size(Z))
            error('p.BC_dir_value should be a (sparse) matrix with the same dimensions as the modelled domain (DEM.Z)');
        end
        switch BC.nbGhost
            case 1
                Z([ BC.BC_indices.topRow  BC.BC_indices.botRow  BC.BC_indices.leftRow  BC.BC_indices.rightRow])=BC.dir_value([ BC.BC_indices.topRow  BC.BC_indices.botRow  BC.BC_indices.leftRow  BC.BC_indices.rightRow]);
            case 2
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2])=BC.dir_value([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2]);
        end
    case 'Neumann'
        %Currently, this is implemented only for a constant, zero gradient at the boundary cells.
        switch BC.nbGhost
            case 1
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow  BC.BC_indices.botRow  BC.BC_indices.topRow  ])=Z([ BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow_2   BC.BC_indices.topRow_2 ]);
            case 2
                Z([ BC.BC_indices.rightRow  BC.BC_indices.leftRow   BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2  BC.BC_indices.botRow  BC.BC_indices.topRow  BC.BC_indices.botRow_2  BC.BC_indices.topRow_2 ])=Z([ BC.BC_indices.rightRow_2  BC.BC_indices.leftRow_2 BC.BC_indices.rightRow_3  BC.BC_indices.leftRow_3  BC.BC_indices.botRow_2   BC.BC_indices.topRow_2  BC.BC_indices.botRow_3  BC.BC_indices.topRow_3 ]);
        end
    otherwise
        error('BC not properly set!')
end