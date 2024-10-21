clc;
clear all;
close all;

% input for the example in the lecture notes
E = 70e3; % E-modulus
BC = [1 2 4 5 6 8; 0 0 0 0 0 0 ]; %boundary conditions; form: [degrees of freedom; applied displacements]
F = [0 0 1500 0 0 0 -500 0 0 0]; % force vector
Con = [1 2; 3 5 ; 4 5 ; 5 1];  %connectivity matrix; form: [element 1; element 2;etc.]
NodePos = [1250 750; 2000 0 ; 750 0 ; 0 0; 500 500]; % node positions; form: [node 1; node 2; etc]
Area = [50 30 30 50]; % vector with the area for each element

% functions here
[u, strain, stress, reactions] = PerseveranceFEMfunctions(E,BC,F,Con,NodePos,Area);

% output: do not change!
for el= 1:size(Con,1)
	fprintf(['the stress in element ', num2str((el)) ' equals ', num2str(stress(el)), ' MPa \n'])
end

for node=1:size(NodePos,1)
	fprintf(['The displacement of node ', num2str(node) , ' in x-direction is ', num2str(u(2*(node)-1)), ' mm \n'])
	fprintf(['The displacement of node ', num2str(node) , ' in y-direction is ', num2str(u(2*(node))), ' mm \n'])
	fprintf(['The reaction force at node ', num2str(node) , ' in x-direction is ', num2str(reactions(2*(node)-1)), ' N \n'])
	fprintf(['The reaction force at node ', num2str(node) , ' in y-direction is ', num2str(reactions(2*(node))), ' N \n']) 
end