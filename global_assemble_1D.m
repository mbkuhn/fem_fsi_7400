function [gradG,resG] = global_assemble_1D(gradG,resG,gradGl,resGl,econn,J,BC)
% Assembles the global stiffness matrix and force vector
% based on the element information supplied.

npe = length(econn);

for nna=1:npe
    A = econn(nna);
    % Boundary conditions
    BCflag = 0;
    [nBC,~] = size(BC);
    for iBC=1:nBC
        % Check node list for BC
        if (ismember(A,BC{iBC,3}))
            % If dirichlet, allow no change to node value
            if (BC{iBC,1}=='d')
                resG(A) = 0;
                gradG(A,:) = 0;
                gradG(A,A) = 1;
                BCflag = 1;
            end
        end
    end
    
    % The rest of the domain
    if (BCflag == 0)
        % Global force vector
        resG(A) = resG(A) + resGl(nna)*J; % in x
        for nnb=1:npe
            B = econn(nnb);
            % Global stiffness matrix
            gradG(A,B) = gradG(A,B) + gradGl(nna,nnb)*J; % x,x
        end
    end
end