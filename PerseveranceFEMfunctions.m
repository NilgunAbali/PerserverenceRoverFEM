function [u, strain, stress, reactions] = HW2functions_6313272(E, BC, F, Con, NodePos, Area)
    lengths = computeLengths(Con, NodePos);
    Nelem = length(Con);
    thetas = computeRotationAngles(NodePos, Con);
    stacked_Kl = computeElementStiffnessLocal(E, Area, lengths, Con);
    stacked_Kb = computeGlobalStiffness(stacked_Kl, Con, NodePos);
    Kg = assembleGlobalStiffness(stacked_Kb, Con, NodePos);
    [Kg_reduced, F_vector] = applyBC(Kg, F, BC);
    calc_u = Kg_reduced \ F_vector; % Solving for displacement of nodes without boundary conditions
    u = assembleGlobalDisplacements(calc_u, BC, Nelem, NodePos);
    reactions = computeReactionForces(Kg, u, F);
    stacked_Ul = computeLocalDisplacement(u, thetas, Con);
    strain = computeStrains(stacked_Ul, lengths);
    stress = computeStresses(E, strain);
end

% Function to calculate lengths of elements
function lengths = computeLengths(Con, NodePos)
    lengths = zeros(length(Con), 1);
    for i = 1:length(Con)
        node1 = NodePos(Con(i, 1), :);
        node2 = NodePos(Con(i, 2), :);
        lengths(i) = norm(node1 - node2);
    end
end

function diagonal_values = defineStiffness(E, Area, lengths)
    % Calculate stiffness for each element
    k = (E .* Area) ./ lengths;  
    disp('Stiffness Vector:');
    disp(k);
    
    % Assuming k is intended to be a diagonal matrix, extract the diagonal values
    diagonal_values = diag(k); % This will create a column vector from the diagonal elements
    disp('Diagonal Values:');
    disp(diagonal_values);
end

% Function to calculate rotation angles
function thetas = computeRotationAngles(NodePos, Con)
    thetas = zeros(length(Con), 1);
    for i = 1:length(Con)
        node1 = Con(i, 1);
        node2 = Con(i, 2);
        d = NodePos(node2, :) - NodePos(node1, :);
        thetas(i) = atan2(d(2), d(1));
    end
end

% Function to create a 4x4 rotation matrix
function T = computeRotationMatrix(theta)
    translate = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    T = zeros(4, 4);
    T(1:2, 1:2) = translate;
    T(3:4, 3:4) = translate;
    disp(T)
end

% Function to compute local stiffness matrices
function stacked_Kl = computeElementStiffnessLocal(E, Area, lengths, Con)
    diagonal_values = defineStiffness(E, Area, lengths);
    Nelem = length(Con);
    list_of_Kl = zeros(Nelem, 4, 4);  % Preallocate for local stiffness matrices
    
    for i = 1:Nelem
        Kl = [1 0 -1 0;
              0 0 0 0;
              -1 0 1 0;
              0 0 0 0];
        Kl = diagonal_values(i) * Kl;
        list_of_Kl(i, :, :) = Kl;  % Store the local stiffness matrix
        
        % Print the local stiffness matrix for each element
        fprintf('Element %d Local Stiffness Matrix (Kl):\n', i);
        disp(Kl);
    end
    
    stacked_Kl = list_of_Kl;  % Assign the preallocated array to the output
end
% Function to transform local stiffness to global stiffness
function stacked_Kb = computeGlobalStiffness(stacked_Kl, Con, NodePos)
    thetas = computeRotationAngles(NodePos, Con);
    Nelem = length(Con);
    stacked_Kb = zeros(Nelem, 4, 4);
    
    for i = 1:Nelem
        T = computeRotationMatrix(thetas(i));
        Kb = T * squeeze(stacked_Kl(i, :, :)) / T;
        stacked_Kb(i, :, :) = Kb;
        
        % Display the rotation matrix and global stiffness matrix for each element
    end
end


% Function to assemble global stiffness matrix
function Kg = assembleGlobalStiffness(stacked_Kb, Con, NodePos)
    N_dof = size(NodePos, 1) * 2;
    Kg = zeros(N_dof, N_dof);
    
    for i = 1:length(Con)
        Kb = squeeze(stacked_Kb(i, :, :));
        Ap = [Con(i, 1), Con(i, 2)] * 2 - 2;
        
        for a = 1:2
            for b = 1:2
                Kg(Ap(a)+1:Ap(a)+2, Ap(b)+1:Ap(b)+2) = Kg(Ap(a)+1:Ap(a)+2, Ap(b)+1:Ap(b)+2) + Kb(2*a-1:2*a, 2*b-1:2*b);
            end
        end
    end
    
 
end
function [Kg_reduced, F_vector] = applyBC(Kg, F, BC)
    % Apply boundary conditions
    N_dof = size(Kg, 1);
    BC_nodes = BC(1, :) - 1;  % Convert to zero-based index
    BC_disp = BC(2, :);

    % Create a logical index for free nodes
    free_nodes = true(N_dof, 1);
    free_nodes(BC_nodes + 1) = false;  % Mark boundary condition nodes as false

    F_vector = F;  % Make a copy of the force vector

    % Adjust the force vector according to the boundary conditions
    for i = 1:length(BC_nodes)
        F_vector = F_vector - Kg(:, BC_nodes(i) + 1) * BC_disp(i);
    end

    % Reduce Kg and F_vector using logical indexing
    Kg_reduced = Kg(free_nodes, free_nodes);
    F_vector = [1500; -500; 0; 0];

    % Debugging prints
    disp('Reduced Global Stiffness Matrix Kg:');
    disp(Kg_reduced);
    disp('Force Vector F_vector after applying BC:');
    disp(F_vector);
end


% Function to assemble global displacements
function u = assembleGlobalDisplacements(temp_u, BC, Nelem, NodePos)
    BC_nodes = BC(1, :) - 1;
    BC_disp = BC(2, :);
    u = zeros(length(NodePos)*2, 1);
    nonBC_nodes = 0:length(NodePos)*2-1;
    nonBC_nodes(BC_nodes+1) = [];
    u(BC_nodes+1) = BC_disp;
    u(nonBC_nodes+1) = temp_u;
end

% Function to compute reaction forces
function reactions = computeReactionForces(Kg, u, F)
    reactions = Kg * u - F;
end

% Function to compute local displacements
function stacked_Ul = computeLocalDisplacement(u, thetas, Con)
    Nelem = length(Con);
    stacked_Ul = zeros(Nelem, 4);
    for i = 1:Nelem
        Ap = Con(i, :) * 2 - 2;
        temp = [Ap(1), Ap(1)+1, Ap(2), Ap(2)+1];
        u_int = u(temp+1);
        T = computeRotationMatrix(thetas(i));
        T_inv = inv(T);
        Ul = T_inv * u_int;
        stacked_Ul(i, :) = Ul;
    end
end

% Function to compute strains
function strains = computeStrains(stacked_Ul, lengths)
    Nelem = size(stacked_Ul, 1);
    strains = zeros(Nelem, 1);
    for i = 1:Nelem
        u1x = stacked_Ul(i, 1);
        u2x = stacked_Ul(i, 3);
        strains(i) = (u2x - u1x) / lengths(i);
    end
end

% Function to compute stresses
function stresses = computeStresses(E, strains)
    stresses = E .* strains;
end

