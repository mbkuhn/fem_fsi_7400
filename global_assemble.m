function [gradG,resG] = global_assemble(gradG,resG,gradGl,resGl,econn,J,BC)
% Assembles the global stiffness matrix and force vector
% based on the element information supplied.

npe = length(econn);
for nna=1:npe
    A = econn(nna);
    
    Rind = 3*A-2:3*A; rind = 3*nna-2:3*nna;
    
    % Boundary conditions
    BCflag = zeros(3,1);
    [nBC,~] = size(BC);
    for iBC=1:nBC
        % Check node list for BC
        if (fastintersect(A,BC{iBC,3}))
            % If dirichlet, allow no change to node value
            if (BC{iBC,1}(1)=='d')
                if (BC{iBC,1}(2)=='x')
                    ind = 3*A-2;
                    BCflag(1) = 1;
                elseif (BC{iBC,1}(2)=='y')
                    ind = 3*A-1;
                    BCflag(2) = 2;
                else % (BC{iBC,1}(2)=='p')
                    ind = 3*A;
                    BCflag(3) = 3;
                end
                resG(ind) = 0;
                gradG(ind,:) = 0;
                gradG(ind,ind) = 1;
            end
        end
    end
    
    % Remove zero values from BCflag
    BCf1 = BCflag(BCflag~=0);
    % According to BCflag, skip certain indices
    Rind(BCf1) = []; rind(BCf1) = [];    
    
    % Global force vector
    resG(Rind) = resG(Rind) + resGl(rind)*J;
    for nnb=1:npe
        B = econn(nnb);
        
        Cind = 3*B-2:3*B; cind = 3*nnb-2:3*nnb;
        % Global stiffness matrix
        gradG(Rind,Cind) = gradG(Rind,Cind) ...
            + gradGl(rind,cind)*J;
    end
    
end