function [u0,udot0] = set_initial_state_1D(mesh,BC)
% Assembles the global stiffness matrix and force vector
% based on the element information supplied.

% Get number of nodes
nnodes = max(max(mesh.conn));

% Initialize size
u0 = zeros(nnodes,1);
udot0 = u0;

% Loop through node list of each dirichlet boundary condition

[nBC,~] = size(BC);
for iBC=1:nBC
    % If dirichlet, set initial state according to BC
    if (BC{iBC,1}=='d')
        for n=1:length(BC{iBC,3})
            ind = BC{iBC,3}(n);
            u0(ind) = BC{iBC,2};
        end
    end
end