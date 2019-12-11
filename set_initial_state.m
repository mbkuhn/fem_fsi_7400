function [u0,udot0,p0,yddot0,ydot0,y0] = set_initial_state(mesh,BC,init_var)
% Assembles the global stiffness matrix and force vector
% based on the element information supplied.

% Get number of nodes
nnodes = mesh.nnode;

% Initialize size
p0 = zeros(nnodes,1);
u0 = zeros(2*nnodes,1);
udot0 = u0;
yddot0 = u0; ydot0 = u0; y0 = u0;

% Initial values, as designated by user
u0(1:2:2*nnodes-1) = init_var(1);
u0(2:2:2*nnodes  ) = init_var(2);
y0(1:2:2*nnodes-1) = init_var(3);
y0(2:2:2*nnodes  ) = init_var(4);

% Loop through node list of each dirichlet boundary condition

[nBC,~] = size(BC);
for iBC=1:nBC
    % If dirichlet, set initial state according to BC
    if (BC{iBC,1}(1)=='d')
        for n=1:length(BC{iBC,3})
            ind = BC{iBC,3}(n);
            if (BC{iBC,1}(2)=='x')
                if (BC{iBC,1}(3)=='f') % fluid
                    u0(2*ind-1) = BC{iBC,2};
                else % (BC{iBC,1}(3)=='s'||BC{iBC,1}(3)=='m') solid or mesh
                    y0(2*ind-1) = BC{iBC,2};
                end
            elseif (BC{iBC,1}(2)=='y')
                if (BC{iBC,1}(3)=='f') % fluid
                    u0(2*ind  ) = BC{iBC,2};
                else % (BC{iBC,1}(3)=='s'||BC{iBC,1}(3)=='m') solid or mesh
                    y0(2*ind  ) = BC{iBC,2};
                end
            else % (BC{iBC,1}(2)=='p')
                p0(  ind  ) = BC{iBC,2};
            end
        end
    end
end

% Zero velocity of solid nodes 
% (so that the initial u doesn't show up in solid)
for iel=1:mesh.nelem
    if (mesh.tag(iel)==1)
        for nna=1:mesh.npe
            ind=mesh.conn(iel,nna);
            u0(2*ind-1:2*ind) = 0.0;
        end
    end
end